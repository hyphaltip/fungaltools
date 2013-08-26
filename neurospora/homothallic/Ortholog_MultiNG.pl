#!/usr/bin/perl -w
use strict;
use Bio::DB::SeqFeature::Store;
use Getopt::Long;
use Env qw(HOME);
# here is cmdline example
# annot_xfer_mercatorblocks.pl -f gene:Broad -t gene:maker \
# -dt Neurospora_pannonica \
# -df Neurospora_crassa -m alignments \
# -gt npan -gf ncra > npan-npan.mercator_orthologs.tab

my ($user,$pass,$dbname_from,@dbname_to,$host);
$host ='localhost';
my $prefix;
my $debug = 0;
my $to_feature = 'gene:maker';
my $from_feature = 'gene:maker';

my ($mercator_dir,$output);
my ($genome_from,@genome_to);

GetOptions(
'v|verbose!' => \$debug,
'u|user:s' => \$user,
'p|pass:s' => \$pass,
'host:s' => \$host,
'df|dbfrom:s' => \$dbname_from,
'dt|dbto:s' => \@dbname_to,

'f|from:s' => \$from_feature,
't|to:s' => \$to_feature, 
'gf|genomefrom:s' => \$genome_from,
'gt|genometo:s' => \@genome_to,

'm|mercator:s' => \$mercator_dir, # alignment dir
'output:s' => \$output,
);

unless( defined $dbname_from ) {
    die("no dbname_from provided\n");
}

unless(  @dbname_to ) { #deleted defined,, but still bugs -related to structure?

    die("no dbname_to provided\n");
}

unless( defined $mercator_dir && -d $mercator_dir ) {
    die("cannot open $mercator_dir, provide with -m or --mercator\n");
}

if( $output && $output ne '-' ) {
    open($output => ">$output" ) || die $!;
} else {
    $output = \*STDOUT;
}

unless( defined $genome_from ) {
    die("must provide a query genome name with -gf or --genomefrom\n");
}

unless( @genome_to ) {
    die("must provide  query genomes name with -gt or --genometo\n");
}

($user,$pass) = &read_cnf($user,$pass) unless $pass && $user;
my $dsn = sprintf('dbi:mysql:database=%s;host=%s',$dbname_from,$host);
my $dbh_from = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::mysql',
-dsn => $dsn,
-user => $user,
-password => $pass,
);


my %handles;

for my $name ( @dbname_to ) {
  $dsn = sprintf('dbi:mysql:database=%s;host=%s',$name,$host);
  my $dbh_to = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::mysql',
					       -dsn => $dsn,
					       -user => $user,
					       -password => $pass,
					      );
  $handles{$name} = $dbh_to;  #$dsn_to
}

my $iter = $dbh_from->get_seq_stream(-type => $from_feature);
my (undef,$fromsrc) = split(/:/,$from_feature);
my (undef,$tosrc) = split(/:/,$to_feature);

print $output join("\t", qw(GENE_FROM MRNA_FROM CHROM_FROM
START_FROM STOP_FROM STRAND_FROM
GENE_TO MRNA_TO START_TO STOP_TO
STRAND_TO SINGLE_MATCH)), "\n";


my %correspondances;
$correspondances{@genome_to} = @dbname_to;  
  
#To use later in replacement of the correspondances WHICH ASSUME THAT dbname_to and genome_to lists contain a listing of species IN THE SAME ORDER: Open a file with correspondances


while( my $gene = $iter->next_seq ) {

    my $name = $gene->name;
    my ($mRNA) = $gene->get_SeqFeatures('mRNA'); # 1st mRNA for now
    if( ! defined $mRNA ) {
      warn("no mRNA for $name\n");
      next;
    }
    my $t_name = $mRNA->name;
    # This program requires the sliceAlignment program from
    # MERCATOR which is part of Colin Dewey's tools
    my $arg = sprintf("sliceAlignment %s %s %s %d %d %s",
$mercator_dir, $genome_from,
$gene->seq_id, $gene->start, $gene->end,
$gene->strand > 0 ? '+' : '-');
    open(my $map => "$arg 2>/dev/null | grep '>' | " ) || die "Cannot open slice with $arg\n";
    my $seen = 0;
    while(<$map>) {
# parsing output from sliceAlignment, all we really want is the
        # matching interval in the other genome (genome_t@ is the prefix)
if( /^>(\S+)\s+([^:]+):(\d+)\-(\d+)([+-])/ ) {
my ($genome,$chrom,$start,$end,$strand) = ($1,$2,$3,$4,$5);

foreach my $genomeIn_tolist (@genome_to) {
  warn( "working on genome: ", $genomeIn_tolist , "\n");
  if( $genome eq $genomeIn_tolist ) {
    $seen ++;
# get that segment in the 'TO' genome
    my  $newname=$correspondances{ $genomeIn_tolist };
    warn( "will search for db: ", $newname , "\n");
    my $segment = $handles{$newname}->segment($chrom,$start,$end);

    if( ! defined $segment ) { 
      warn("cannot find segment in the DB for $chrom:$start..$end ($genome)\n");
      next;
    }

# extract the gene(s) in this interval
    my @genes = $segment->features(-type => $to_feature);
    if( @genes ) {
      for my $g ( @genes ) {
	my ($to_mRNA) = $g->get_SeqFeatures('mRNA'); # +++
# print out the genes that fall in syntenic interval
	  print $output join("\t",
			     $name, $t_name,$gene->seq_id,
			     $gene->start,$gene->end,
			     $gene->strand,

			     $g->name,$to_mRNA->name,
			     $g->seq_id,$g->start,$g->end,
			     $g->strand,
			     @genes == 1 ? 'yes' : 'no'),"\n";
	}
      } else {
	print $output join("\t",
			   $name, $t_name,$gene->seq_id,
			   $gene->start,$gene->end,$gene->strand,
			   '','','','','','NO_GENES_IN_INTERVAL'),"\n";
      }

}
}
    }
    unless( $seen ) { # To edit!!!
print $output join("\t",
$name, $t_name,$gene->seq_id,
$gene->start,$gene->end,$gene->strand,
'','','','','','NO_SYNTENIC_ALIGNMENT'),"\n";
    }
    last if $debug;
}
  }

# read the .my.cnf file
sub read_cnf {
    my ($user,$pass) = @_;
    if( -f "$HOME/.my.cnf") {
        open(IN,"$HOME/.my.cnf");
        while(<IN>) {
            if(/user(name)?\s*=\s*(\S+)/ ) {
                $user = $2;
            } elsif( /pass(word)\s*=\s*(\S+)/ ) {
                $pass = $2;
            }
        }
        close(IN);
    }
    return ($user,$pass);
}
