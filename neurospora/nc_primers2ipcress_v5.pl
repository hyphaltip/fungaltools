#!/usr/bin/perl -w
use strict;
# Author: Jason Stajich jason.stajich@ucr.edu
use strict;

# based on this file,converted to tab delimited
# from the ver5 file provided by Carol
# neurospora_crassa_or74a_finished_10_supercontigs_new_primers.txt

my $in = 'neurospora_crassa_or74a_finished_10_supercontigs_new_primers.txt';

open(my $left => ">primer5.ver5.dat");
open(my $right => ">primer3.ver5.dat");

# 5primer codes are
my $leftremove5   = 'GTAACGCCAGGGTTTTCCCAGTCACGACG';
#my $rightremove5 = 'ATCCACTTAACGTTACTGAAATCTCCAAC';
my $rightremove5 = 'ATCCACTTAACGTTACTGAAATCT';

# 3prime codes are
my $leftremove3  = 'CTCCTTCAATATCATCTTCTGTCTCCGAC';
my $rightremove3 = 'GCGGATAACAATTTCACACAGGAAACAGC';

open(my $fh => $in) || die $!;
while(<$fh>) {
    next if /^locus/ || /^\s+/;
    last if /^Summary/;
    my @row = split;
    my ($gene,$strand,$start,$end, 
	$status5,
	$flanklen5,$size5,
	$gap5,
	$forward5start,
	$primer5f,
	$reverse5start,
	$primer5r,
	$status3,
	$flanklen3,$size3,
	$gap3,
	$forward3start,
	$primer3f,
	$reverse3start,
	$primer3r,
	) = split;
    
    if( $status5 eq 'failed' || $status3 eq 'failed') {
	print("skipping $gene was it failed\n");
	next;
   }

    if( $primer5f ) {
	$primer5f = substr($primer5f,length($leftremove5));
	$primer5r = substr($primer5r,length($rightremove5));
	print $left join("\t",$gene.".5prime",
			 $primer5f, $primer5r, 
			 $size5 - 300, $size5 + 300),"\n";
    }
    if( $primer3f ) {
	$primer3f = substr($primer3f,length($leftremove3));
	$primer3r = substr($primer3r,length($rightremove3));
	print $right join("\t",$gene.".3prime",
			  $primer3f, $primer3r, $size3 - 300,
			  $size3 + 300),"\n";
    }    
    
}
