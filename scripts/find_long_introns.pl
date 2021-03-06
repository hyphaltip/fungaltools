#!/usr/bin/perl -w

=head1 NAME

find_long_introns

=head1 DESCRIPTION

Fungi typically have short introns < 500 bp so long introns may indicate 
misannotation

=head1 AUTHOR

Jason Stajich jason.stajich[at]ucr.edu

=cut

use strict;
use warnings;

use Bio::DB::SeqFeature::Store;
use Getopt::Long;
use Env qw(HOME);
my ($user,$pass,$dbname,$host);
$host ='localhost';
my $prefix;
my $debug = 0;
my $src = 'gene:NC10_CALLGENES_FINAL_5';
my $output;
my $min_intron_size_reported = 300;
GetOptions(
           'v|verbose!'  => \$debug,
           'u|user:s'    => \$user,
           'p|pass:s'    => \$pass,
           'host:s'      => \$host,
           'db|dbname:s' => \$dbname,

           's|src:s'     => \$src,
           'o|output:s'  => \$output,
	   'm|min:i'     => \$min_intron_size_reported,
           );

unless(  defined $dbname ) {
    die("no dbname provided\n");
}

($user,$pass) = &read_cnf($user,$pass) unless $pass && $user;
my $dsn = sprintf('dbi:mysql:database=%s;host=%s',$dbname,$host);
my $dbh = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::mysql',
                                          -dsn     => $dsn,
                                          -user    => $user,
                                          -password => $pass,
                                          );
my $ofh;
if( $output && $output ne '-' ) { 
    open($ofh => ">$output" ) || die $!;
} else {
    $ofh = \*STDOUT;
}


my $iter = $dbh->get_seq_stream(-type => $src);
my $count = 0;
my $gene_count = 0;
my @introns;
while( my $gene = $iter->next_seq ) {
    $gene_count++;
    my $gene_name = $gene->name;
    my $sourcetag = $gene->source;
    my @mRNA = $gene->get_SeqFeatures("mRNA:$sourcetag");
    if( ! @mRNA ) {
        warn("no mRNA for $gene_name\n");
    }
    for my $mRNA ( @mRNA ) { # 1st mRNA for now
        my $last_exon;
        my $i = 1;
        my $mRNA_name = $mRNA->name || 'mRNA-'.$mRNA->load_id;
        for my $exon ( sort { $a->start * $a->strand <=> 
                                  $b->start * $b->strand } 
                       $mRNA->get_SeqFeatures("exon:$sourcetag") ) {
            if( $last_exon ) {
                my ($start,$end) = ( $last_exon->end+1,
                                     $exon->start - 1);
                if( $exon->strand < 0 ) {
                    ($start,$end) = ( $last_exon->start-1,
                                      $exon->end + 1);
                }
		push @introns, [$gene->seq_id,
                                $start,$end,
                                $exon->strand,
				abs($end-$start),
                                sprintf('ID=%s.i%d;Gene=%s',
                                        $mRNA_name,
                                        $i++,
                                        $gene_name)];
	      }
            $last_exon = $exon;
        }
        last;
    }
    last if $debug && $count++ > 10;
}
warn("gene count $gene_count\n");

for my $intron ( sort { $b->[4] <=> $a->[4] }
		 @introns
	       ) {
  last if $intron->[4] < $min_intron_size_reported;
  print join("\t", @$intron),"\n";
}

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



