#!/usr/bin/perl -w
use strict;

use DBI;
use Getopt::Long;

my %uncomp = ('bz2' => 'bzcat',
	      'gz'  => 'gzcat');

my $interval = 100_000;

my $db = 'pairs.db';
my $force = 0;
my $input;
my $auto = 0;
my $debug = 0;
my $cutoff_id = 20;
my $cutoff_pmatch = 30;

my $outgroup_filter = 'Hsap';

my (@ignore);
GetOptions(
	   'd|db:s' => \$db,
	   'force!' => \$force,
	   'i|in:s' => \$input,
	   'outgroup:s' => \$outgroup_filter,
	   'commitsize:i' => $interval,
	   'v|debug!' => \$debug,
	   'ignore=s' => \@ignore,
);

my %ignore = map { $_ => 1 } @ignore;

$input= shift @ARGV unless defined $input;

my $exists = (-f $db && ! -z $db);
my $dbh = DBI->connect("dbi:SQLite:dbname=$db","","", 
		       {AutoCommit => $auto});
#  $dbh->do("PRAGMA foreign_keys=ON");

  my $i = 0;

if ( $force || ! $exists ) {
  unlink($db);
  $dbh = DBI->connect("dbi:SQLite:dbname=$db","","", 
		      {AutoCommit => $auto});
  #  $dbh->do("PRAGMA foreign_keys=ON");
  $dbh->do(&create_tables());
  my $sth = $dbh->prepare("INSERT INTO SimilarSequences VALUES (?,?,?,?,?,?,?,?)");
  my $fh;

  if ( $input =~ /\.(bz2|gz|Z)/ ) {
    open($fh => "$uncomp{$1} $input |") || die $!;
  } else {
    open($fh => $input) || die $!;
  }

  while (<$fh>) {
    $sth->execute(split);
    $dbh->commit if $i++ % $interval == 0;
  }
  $dbh->commit;
}

my %species;

my $qsth = $dbh->prepare("SELECT query_taxon_id, COUNT(DISTINCT query_id) FROM SimilarSequences GROUP BY query_taxon_id");

$qsth->execute();
my ($qname,$ct) ;
$qsth->bind_columns(\$qname, \$ct);
while ($qsth->fetch ) {
  $species{$qname} = $ct;
}
$qsth->finish;

$qsth = $dbh->prepare("SELECT DISTINCT query_id FROM SimilarSequences WHERE query_taxon_id = ?");
my $qc_sth = $dbh->prepare(<<EOF
SELECT subject_taxon_id, subject_id, evalue_exp
FROM SimilarSequences WHERE query_id = ? AND percent_identity >= ? AND percent_match >= ? AND subject_id != ? ORDER BY subject_taxon_id
EOF
);

$i = 0;
my %patterns;
my %gene;
my %genelist;

for my $sp ( keys %species ) {
  next if $ignore{$sp};
  $qsth->execute($sp);
  my ($subject_taxid,$subject, $evalue, $pid, $pmatch);
  $qsth->bind_columns(\$qname);
  while ($qsth->fetch ) {
    $qc_sth->execute($qname, $cutoff_id, $cutoff_pmatch,$qname);
    $qc_sth->bind_columns(\$subject_taxid,\$subject,\$evalue );
    my %snames = ($sp => [0,$qname] );
    while($qc_sth->fetch ) {
      # skip if this is an ignore taxa
      next if $subject eq $qname;
      next if $ignore{$subject_taxid};
      if ( ! exists $snames{$subject_taxid} ) {
	$snames{$subject_taxid} = [$evalue,$subject];
#	warn("$subject_taxid = $evalue - $qname\n");
      } elsif ( $snames{$subject_taxid}->[0] != 0 &&
		$snames{$subject_taxid}->[0] > $evalue ) {
#	warn("replacing $snames{$subject_taxid} with $evalue for $subject_taxid and $qname\n");
	$snames{$subject_taxid} = [$evalue,$subject];
      }
    }
    # filter by outgroup (human)
    # to remove hits when the outgroup has a better Evalue
    if ( exists $snames{$outgroup_filter} ) {
      my @to_remove;
      while ( my ($taxid,$tax_evalue) = each %snames ) {
	next if $outgroup_filter eq $taxid || $tax_evalue->[0] == 0;
	if ( ($tax_evalue->[0] > $snames{$outgroup_filter}->[0]) &&
	     abs($tax_evalue->[0] / $snames{$outgroup_filter}->[0])  < 0.20 ) {
	  warn("removing ",join("-",@{$tax_evalue}),
	       " for (outgroup is ",
	       join('-',@{$snames{$outgroup_filter}}),,") for $qname\n");
	  push @to_remove, $taxid;
	}
      }
      # remove taxa where the evalue is worse than the outgroup one
      for my $r ( @to_remove ) {
	delete $snames{$r};
      }
#      if ( @to_remove ) {
#	warn("removing @to_remove for $qname\n");
#	exit;
#      }
    }
    warn join("\t", $qname, sort keys %snames), "\n" if $debug;
    my $pat = join(",", sort keys %snames);
    push @{$patterns{$sp}->{$pat}}, $qname;
    $gene{$sp}->{$qname}->{pat} = $pat;
    $gene{$sp}->{$qname}->{hits} = \%snames;
    $qc_sth->finish;
    last if $debug && $i++ > 100;
  }
  $qsth->finish;
  last if $debug;
}
$dbh->commit;
$dbh->disconnect;

for my $sp ( keys %patterns ) {
  open( my $ofh => ">$sp.patterns.tab") || die $!;
  open( my $ofhg => ">$sp.gene_patterns.tab") || die $!;
  open( my $ofhgn => ">$sp.gene_clusters.tab") || die $!;

  # this mapping will mean that $p has in it two values in the arrayref pattern and count
  for my $p ( sort { $b->[1] <=> $a->[1] }  # sort by number of members in each pattern
	      map { [ $_, scalar @{$patterns{$sp}->{$_}} ] }
	      keys %{$patterns{$sp}} ) {
    my ($pat,$pcount) = @$p;
    print $ofh join("\t", $pcount, $pat), "\n";
    for my $g ( @{$patterns{$sp}->{$pat}} ) {
      print $ofhg join("\t", $g, $pat), "\n";
    }
  }

  for my $gn ( keys %{$gene{$sp}} ) {
    my $hits = $gene{$sp}->{$gn}->{hits};
    print $ofhgn join("\t", $gn, map { sprintf("%s=%s",(reverse @{$hits->{$_}})) } keys %{$hits}), "\n";
  }
  close($ofh);
}

sub create_tables {

my $sql = <<EOF
CREATE TABLE IF NOT EXISTS SimilarSequences  (
 QUERY_ID                 VARCHAR(60),
 SUBJECT_ID               VARCHAR(60),
 QUERY_TAXON_ID           VARCHAR(40),
 SUBJECT_TAXON_ID         VARCHAR(40),
 EVALUE_MANT              FLOAT,
 EVALUE_EXP               INT,
 PERCENT_IDENTITY         FLOAT,
 PERCENT_MATCH            FLOAT
);


CREATE UNIQUE INDEX IF NOT EXISTS ss_qtaxexp_ix
ON SimilarSequences(query_id, subject_taxon_id,
evalue_exp, evalue_mant,
query_taxon_id, subject_id);

CREATE UNIQUE INDEX IF NOT EXISTS ss_seqs_ix
ON SimilarSequences(query_id, subject_id,
evalue_exp, evalue_mant, percent_match);

CREATE INDEX IF NOT EXISTS qtaxon_id
ON SimilarSequences(query_taxon_id);

CREATE INDEX IF NOT EXISTS staxon_id
ON SimilarSequences(subject_taxon_id);

CREATE INDEX IF NOT EXISTS q_id
ON SimilarSequences(query_id);

CREATE INDEX IF NOT EXISTS s_id
ON SimilarSequences(subject_id);


EOF
;

$sql;

}
