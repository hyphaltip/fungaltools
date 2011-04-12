#!/usr/bin/perl -w
# Author: Jason Stajich jason.stajich@ucr.edu
use strict;

# based on this file,converted to tab delimited
# http://www.dartmouth.edu/~neurosporagenome/Projects_files/Project%201/primers_for_1kb_flanks.xls

open(my $left => ">primer5.dat");
open(my $right => ">primer3.dat");

# 5primer codes are
my $leftremove5   = 'GTAACGCCAGGGTTTTCCCAGTCACGACG';
my $rightremove5 = 'ATCCACTTAACGTTACTGAAATCTCCAAC';

# 3prime codes are
my $leftremove3  = 'CTCCTTCAATATCATCTTCTGTCTCCGAC';
my $rightremove3 = 'GCGGATAACAATTTCACACAGGAAACAGC';

while(<>) {
    next if /^revised/ || /^GENE/;
    chomp;
    my ($gene,$primer5f, $primer5r, $primer3f, $primer3r,$gap5,$gap3,$size5,$size3) = split(/\t/,$_);
    if( $primer5f ) {
	$primer5f = substr($primer5f,length($leftremove5));
	$primer5r = substr($primer5r,length($rightremove5));
	print $left join("\t",$gene.".5prime",
			 $primer5f, $primer5r, $size5 - 300, $size5 + 300),"\n";
    }
    if( $primer3f ) {
	$primer3f = substr($primer3f,length($leftremove3));
	$primer3r = substr($primer3r,length($rightremove3));
	print $right join("\t",$gene.".3prime",
			  $primer3f, $primer3r, $size3 - 300, 
			  $size3 + 300),"\n";
    }
}
