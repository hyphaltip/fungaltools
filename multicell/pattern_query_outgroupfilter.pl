#!/usr/bin/perl -w
use strict;
use DBI;
use File::Spec;
use Getopt::Long;
use Env qw(USER);

my %uncomp = ('bz2' => 'bzcat',
	      'gz'  => 'gzcat');

my $interval = 10_000_000;

my $db = 'pairs.db';
my $tmpdir = '/scratch'; # probably should be /tmp on other systems but we use /scratch for fast SSD drives
my $make_in_tmp = 1;
my $force = 0;
my $input;
my $auto = 0;
my $debug = 0;
my $runonce = 0;
my $cutoff_id = 20;
my $cutoff_pmatch = 30;
my $evalue_outgroup_ratio = 0.20;
my $evalue_cutoff = 1e-4;
my $generate_stats;
my $outgroup_filter = 'HSAP';
my $build_only;
my $only_taxon;
my $odir = ".";
my (@ignore);
GetOptions(
	   'd|db:s' => \$db,
	   'force!' => \$force,
	   'i|in:s' => \$input,
           'o|outdir|dir:s' => \$odir,
	   'outgroup:s' => \$outgroup_filter,
	   'commitsize:i' => $interval,
	   'v|debug!' => \$debug,
	   'runonce!' => \$runonce,
	   'or|outratio:s' => \$evalue_outgroup_ratio,
	   'identity:s' => \$cutoff_id,
	   'pmatch:s'   => \$cutoff_pmatch,
	   'showstats!' => \$generate_stats,
	   'e|evalue:s' => \$evalue_cutoff,
	   'ignore=s' => \@ignore,
	   'only_taxon:s' => \$only_taxon,
	   'buildonly!' => \$build_only,
  	   'tmpdir:s'   => \$tmpdir,
);

mkdir($odir) unless -d $odir;
my %ignore = map { $_ => 1 } @ignore;

$input= shift @ARGV unless defined $input;

my $exists = (-f $db && ! -z $db);


  my $i = 0;

if ( $force || ! $exists ) {
  unlink($db);
  my (undef,undef,$dbfile) = File::Spec->splitpath($db);
  my $tmpdb = $db;
  if( $make_in_tmp) {
   $tmpdb = File::Spec->catfile($tmpdir,"$USER.$$.$db");
  }
  my $dbh = DBI->connect("dbi:SQLite:dbname=$tmpdb","","",
			 {AutoCommit => 0, RaiseError => 1});
  #  $dbh->do("PRAGMA foreign_keys=ON");
  if( &create_tables($dbh) != 0 ) {
      warn("Error in creating tables. Aborting.\n");
      $dbh->rollback;
      exit;
  }
  my $sth = $dbh->prepare("INSERT INTO SimilarSequences VALUES (?,?,?,?,?,?,?,?)");
  my $fh;

  if ( $input =~ /\.(bz2|gz|Z)/ ) {
    open($fh => "$uncomp{$1} $input |") || die $!;
  } else {
    open($fh => $input) || die "cannot open $input: $!";
  }
  eval {
      while (<$fh>) {
	  my @row = split;
	  $sth->execute(@row);
	  $dbh->commit if $i++ % $interval == 0;
      }
  };
  if( $@ ) {
      warn("$@");
      $dbh->rollback;
      exit;
  }

  $dbh->commit;  
  $sth->finish;
  if( &create_indexes($dbh) != 0 ) {
      warn("error in creating indexes, aborting.\n");
      exit;
  }
  $dbh->disconnect;
  if( $make_in_tmp) {
   `rsync -a --progress $tmpdb $db`;
   unlink($tmpdb);
  }
} else { 
    warn("skip no need to build the db file already exists\n");
}
#  $dbh->do("PRAGMA foreign_keys=ON");

if ( $build_only ) {
 warn"stopping early since build DB only was called\n";
 exit;
}
my $dbh = DBI->connect("dbi:SQLite:dbname=$db","","", 
		       {AutoCommit => 1, RaiseError => 1});

my %species;

my $qsth = $dbh->prepare("SELECT query_taxon_id, COUNT(DISTINCT query_id) FROM SimilarSequences GROUP BY query_taxon_id");

$qsth->execute();
my ($qname,$ct) ;
$qsth->bind_columns(\$qname, \$ct);
while ($qsth->fetch ) {
  $species{$qname} = $ct;
}
$qsth->finish;

my $qsth_self = $dbh->prepare("SELECT evalue_mant, evalue_exp FROM SimilarSequences 
WHERE query_id = ? AND subject_id = ? order by evalue_mant DESC");

$qsth = $dbh->prepare("SELECT DISTINCT query_id FROM SimilarSequences WHERE query_taxon_id = ?");
my $qc_sth = $dbh->prepare(
'SELECT subject_taxon_id, subject_id, evalue_mant, evalue_exp
FROM SimilarSequences WHERE query_id = ? AND percent_identity >= ? AND percent_match >= ? AND subject_id != ?'
);

# ORDER BY subject_taxon_id, subject_id'

$i = 0;
for my $sp ( keys %species ) {
  next if $ignore{$sp};# || $sp eq $outgroup_filter;
  next  if( $only_taxon && $only_taxon ne $sp );
  my $ofile = File::Spec->catfile($odir,"$sp.aln_stats.csv");
  next if -f $ofile && ! -z $ofile;
  my (%patterns,%gene, %genelist);

  open(my $ofh => ">$ofile") || die $!;
  print $ofh join("\t", qw(QUERY SUBJECT EVALUE PERCENT_ID PERCENT_MATCH BRH
			   R_SUBJECT R_QUERY R_EVALUE R_PERCENT_ID R_PERCENT_MATCH)), "\n";

  $qsth->execute($sp);
  my ($subject_taxid,$subject, $e_mant,$e_exp);
  $qsth->bind_columns(\$qname);
  my @qnames;
  while ($qsth->fetch ) {
    $qsth_self->execute($qname,$qname);
    my ($test_m, $test_e);
    $qsth_self->bind_columns(\$test_m, \$test_e);
    unless ( $qsth_self->fetch) {
      warn("error in query for self on $qname vs $qname - no results\n");
      $test_m = 0;
      $test_e = 0;
    }
    my $self_evalue = sprintf("%se%d",$test_m,$test_e);
    $qsth_self->finish;

    my %snames = ($sp => [$self_evalue,$qname] );
    $qc_sth->execute($qname, $cutoff_id, $cutoff_pmatch,$qname);
    $qc_sth->bind_columns(\$subject_taxid,\$subject,\$e_mant, \$e_exp );

    while ($qc_sth->fetch ) {
      next if $subject eq $qname;
      next if $ignore{$subject_taxid};
      my $evalue;
      if ( ! defined $e_mant || ! defined $e_exp ) {
	warn("$qname evalue not found\n");
	$evalue = 0;
      } elsif ( $e_exp == 0 && $e_mant == 0) {
	$evalue = 0;
      } else {
	$evalue = sprintf("%se%d",$e_mant,$e_exp);
      }

      # skip if this is an ignore taxa
      if ( ! exists $snames{$subject_taxid} ) {
	$snames{$subject_taxid} = [$evalue,$subject];
	#	warn("$subject_taxid = $evalue - $qname\n");
      } elsif ( $snames{$subject_taxid}->[0] != 0 &&
		$snames{$subject_taxid}->[0] > $evalue ) {
	  # only keep the BEST hit
	  #	warn("replacing $snames{$subject_taxid} with $evalue for $subject_taxid and $qname\n");
	$snames{$subject_taxid} = [$evalue,$subject];
      }
    }

    # filter by outgroup (human)
    # to remove hits when the outgroup has a better Evalue
    my @to_remove;
    while ( my ($taxid,$tax_evalue) = each %snames ) {
	if( $tax_evalue->[0] > $evalue_cutoff ) {
# if the currently examined evalue is bigger than our current e-value cutoff (e-10) remove this hit
	    warn("removing $taxid as the evalue is ",$tax_evalue->[0], 
		 " while cutoff is ", $evalue_cutoff,"\n");
	    push @to_remove, $taxid;
	    next;
	}
# skip when evalue is 0 or (where sequence is really-really conserved) no one can win in these
# comparisons
	next if( $tax_evalue->[0] == 0 || ! exists $snames{$outgroup_filter} ||
		 ! defined $snames{$outgroup_filter}->[0] ||
		 $snames{$outgroup_filter}->[0] == 0);
	
	if( ($tax_evalue->[0] > $snames{$outgroup_filter}->[0]) &&
	    abs($tax_evalue->[0] / $snames{$outgroup_filter}->[0])  < $evalue_outgroup_ratio ) {
	    # if the evalue is worse to this hit than to the outgroup then
	    if( uc($sp) eq $sp) {
	     warn("removing ",join("-",@{$tax_evalue}),
		 " for (outgroup is ",
		 join('-',@{$snames{$outgroup_filter}}),") for $qname\n") if $debug;
	      push @to_remove, $taxid;
	    } else {
		warn("not removing $taxid for $sp hit because not a fungus\n");
	    }
	}
    }
    # remove taxa where the evalue is worse than the outgroup one or evalue is too large
    for my $r ( @to_remove ) {
	delete $snames{$r};
    }
     
    warn join("\t", $qname, join(",",sort keys %snames)), "\n" if $debug;
    $gene{$sp}->{$qname}->{hits} = \%snames;
    $qc_sth->finish;
#    last if $debug && $i++ > 100;
  }
  $qsth->finish;
  
# generate BRH stats and percent ID stats
  my $sth = $dbh->prepare("SELECT query_id, subject_id, evalue_mant, evalue_exp, percent_identity, percent_match FROM SimilarSequences WHERE query_id = ? and subject_taxon_id = ? ORDER BY evalue_exp, evalue_mant DESC");
  
  for my $gn ( keys %{$gene{$sp}} ) {
    my $hits = $gene{$sp}->{$gn}->{hits};
    my @seensubjtax;
    for my $subjtax ( keys %$hits ) {
      next if $subjtax eq $sp;
      my ($evalue,$sid) = @{$hits->{$subjtax}};
      $sth->execute( $gn, $subjtax); # reverse subject/query
      my ($f_qid, $f_subjid, $f_mant, $f_exp, $f_pid, $f_pmatch);
      while (my $row = $sth->fetchrow_arrayref) {
	if ( $f_qid ) {
	  if ( $row->[3] < $f_exp) {
	    ($f_qid, $f_subjid, $f_mant, $f_exp, $f_pid, $f_pmatch) = @{$row};
	  } else {
	    warn("skipping worse hit for\n\t",join("\t",@$row), " as compared to the better value: \n\t",
		 join("\t",($f_qid, $f_subjid, $f_mant, $f_exp, $f_pid, $f_pmatch)),"\n") if $debug;
	  }
	} else {
	  ($f_qid, $f_subjid, $f_mant, $f_exp, $f_pid, $f_pmatch) = @{$row};
	}
      }
      $sth->finish;
      $sth->execute( $sid,$sp); # reverse subject/query
      my ($r_qid, $r_subjid, $r_mant, $r_exp, $r_pid, $r_pmatch);
      while (my $row = $sth->fetchrow_arrayref) {
	if ( $r_qid ) {
	  if ( $row->[3] < $r_exp) {
	    warn("replacing $r_subjid with ",$row->[1],"\n");
	    ($r_qid, $r_subjid, $r_mant, $r_exp, $r_pid, $r_pmatch) = @{$row};
	  } else {
	    warn("skipping worse hit for\n\t",join("\t",@$row), " as compared to the better value: \n\t",
		 join("\t",($r_qid, $r_subjid, $r_mant, $r_exp, $r_pid, $r_pmatch)),"\n") if $debug;
	  }
	} else {
	  ($r_qid, $r_subjid, $r_mant, $r_exp, $r_pid, $r_pmatch) = @{$row};
	}
      }
      $sth->finish;
      if ( ! $r_qid ) {
	# skip this one the unidirectional search is hitting a false positive it seems
	next;
      }

      my $brh = $f_qid eq $r_subjid ? 'YES' : 'NO';
      print $ofh join("\t", $f_qid, $f_subjid, sprintf("%se%d",$f_mant,$f_exp),
		      $f_pid, $f_pmatch,
		      $brh,
		      $r_qid, $r_subjid, sprintf("%se%d",$r_mant,$r_exp),
		      $r_pid, $r_pmatch), "\n";

      if ( $brh ) {
#	$brhs{$f_qid}->{$f_subjid}++;
#	$brhs{$f_subjid}->{$f_qid}++;
	push @seensubjtax, $subjtax; # require BRH for pattern
      }
    }
    my $pat = join(",", sort @seensubjtax);
    $gene{$sp}->{$gn}->{pats} = $pat;
    push @{$patterns{$sp}->{$pat}}, $gn;
  }
  close($ofh);
  open( $ofh => ">$odir/$sp.patterns.tab") || die $!;
  open( my $ofhg => ">$odir/$sp.gene_patterns.tab") || die $!;
  open( my $ofhgn => ">$odir/$sp.gene_clusters.tab") || die $!;

  # this mapping will mean that $p has in it two values in the arrayref pattern and count
  for my $p ( sort { $b->[1] <=> $a->[1] } # sort by number of members in each pattern
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

  last if $debug && $runonce;
}

$dbh->disconnect;
 
sub create_tables {
    my $dbh = shift;
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

EOF
    ;
    eval { 
	$dbh->do($sql);
	$dbh->commit;
    };
    if( $@ ) {
	warn("Error $@");
	return -1;
    }
    return 0;
}

sub create_indexes {
    my $dbh = shift;

     my @stmts = 
	 ('CREATE UNIQUE INDEX IF NOT EXISTS ss_qtaxexp_ix
ON SimilarSequences(query_id, subject_taxon_id,
evalue_exp, evalue_mant,query_taxon_id, subject_id)',
'CREATE UNIQUE INDEX IF NOT EXISTS ss_seqs_ix
ON SimilarSequences(query_id, subject_id,
evalue_exp, evalue_mant, percent_match)',
'CREATE INDEX IF NOT EXISTS qtaxon_id
ON SimilarSequences(query_taxon_id)',
'CREATE INDEX IF NOT EXISTS staxon_id
ON SimilarSequences(subject_taxon_id)',
'CREATE INDEX IF NOT EXISTS q_id
ON SimilarSequences(query_id)',
'CREATE INDEX IF NOT EXISTS s_id
ON SimilarSequences(subject_id)');

   for my $sql ( @stmts ) {
       eval {
	   my $sth = $dbh->prepare($sql);
	   $sth->execute;
	   $dbh->commit;
       };
       if ( $@ ) {
	   warn($@);
	   $dbh->rollback;
	   return -1;
       }
    }
    return 0;
}


END { 
    $dbh->disconnect if defined $dbh;
}
