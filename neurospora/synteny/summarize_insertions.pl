#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Long;
# summarize the RIP and block length of Reference-only insertions which
# are probably TE relics
my $debug = 0;
GetOptions('v|verbose!' => \$debug);
my $dir = shift;

opendir(DIR, $dir) || die "$dir: $!";
print join("\t", qw(BLOCKNUM CHROM START END STRAND LENGTH RIP_Product RIP_Substrate RIP_CpG)),"\n";
for my $f ( sort { $a->[0] <=> $b->[0] } 
	    map { [/blk_(\d+)\.fasaln/,$_] } 
	    grep { /\.fasaln/ } readdir(DIR) ) {
    my ($stem,$file) = @$f;
    my $seq = Bio::SeqIO->new(-format => 'fasta', -file => "$dir/$file")->next_seq;
    my $len = $seq->length;
    my ($chrom,$start,$end,$strand);
    if( $seq->description =~ /(\S+):(\d+)-(\d+)([+-])/ ){
	($chrom,$start,$end,$strand) = ($1,$2,$3,$4);
    } else {
	warn("cannot parse seq desc ", $seq->description, "\n");
    }	
    print join("\t", $stem, $chrom,$start,$end,$strand, $len, (map { sprintf("%.2f",$_) } rip_index($seq->seq))),"\n";
    last if $debug;
}

use constant DOUBLET_LENGTH => 2;
sub rip_index {
    my $str = shift;
    my $len = length($str);
    my %freq;
    for( my $k = 0; $k < $len - DOUBLET_LENGTH; $k++) {
	my $w = uc(substr($str,$k,DOUBLET_LENGTH));
	$freq{$w}++;
    }
    my $substrate  = ($freq{AC} + $freq{GT}) ?($freq{CA} + $freq{TG})/($freq{AC} + $freq{GT}) : -1;
    return (($freq{AT}) ? $freq{TA}/$freq{AT} : -1, # TpA
	    $substrate, # substrate
	    $freq{CG} ? ($freq{GC} / $freq{CG}) : -1, # GpC
	    );
}
