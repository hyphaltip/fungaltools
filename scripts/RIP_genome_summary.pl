#!/usr/bin/perl -w

# This expects a bed file from RIP_index_calculation.pl

use strict;
use Getopt::Long;
my $cutoff = 0;
my $index= 'CRI';
GetOptions('c|cutoff:s' => \$cutoff,
	   'index:s'    => \$index,
	   );

for my $file ( @ARGV ) {
    open(my $fh => $file) || die "can't open $file: $!";
    
    my $total_windows = 0;
    my $rip_windows = 0;
    while(<$fh>) {
	my ($chrom,$start,$stop,$score) = split(/\s+/,$_);	
	my @row = split;
	if ( $index eq 'CRI' ) {
	  $rip_windows++ if($score > $cutoff );
	} else {
	  $rip_windows++ if($score <= $cutoff );
	}
	$total_windows++;
    }
    printf "%s\t%d windows out of %d total -- %.2f%% RIPed\n",
    $file,$rip_windows, $total_windows, 
    100 * $rip_windows /$total_windows;
}
