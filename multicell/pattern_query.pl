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

my (@ignore);
GetOptions(
	   'd|db:s' => \$db,
	   'force!' => \$force,
	   'i|in:s' => \$input,
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
SELECT DISTINCT subject_taxon_id 
FROM SimilarSequences WHERE query_id = ? AND percent_identity >= ? AND percent_match >= ? ORDER BY subject_taxon_id
EOF
);

$i = 0;
my %patterns;
my %gene;

for my $sp ( keys %species ) {
  next if $ignore{$sp};
  $qsth->execute($sp);
  my ($subject_taxid,$evalue, $pid, $pmatch);
  $qsth->bind_columns(\$qname);
  while ($qsth->fetch ) {
    $qc_sth->execute($qname, $cutoff_id, $cutoff_pmatch);
    $qc_sth->bind_columns(\$subject_taxid );
    my @snames;
    while($qc_sth->fetch ) {
      next if $ignore{$subject_taxid};
      push @snames, $subject_taxid;
    }
    warn join("\t", $qname, @snames), "\n" if $debug;
    my $pat = join(",", sort @snames);
    push @{$patterns{$sp}->{$pat}}, $qname;
    $gene{$sp}->{$qname} = $pat;
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
