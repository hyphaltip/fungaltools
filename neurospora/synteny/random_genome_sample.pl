#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;

my ($datfile,$genome) = @ARGV;

open(my $fh => $datfile) || die $!;
my $dbh = Bio::DB::Fasta->new($genome);
my $header =<$fh>;
print $header;
while(<$fh>) {
    my ($block,$chrom,$start,$end,$strand,$length) = split;
    
    my $chromlen = $dbh->length($chrom);
    my $randstart = int rand($chromlen - $length);

    my $randseq = $dbh->seq($chrom,$randstart+1 => $randstart + $length);
    print join("\t", $block, $chrom, $randstart+1, $randstart + $length,
	       '+', $length, map { sprintf("%.2f",$_) } &rip_index($randseq)),"\n";
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
