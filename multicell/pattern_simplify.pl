#!/usr/bin/perl -w
use strict;
my %groups = ( 'SCOMP' => 'Taphrinamycotina',
	       'SOCT'  => 'Taphrinamycotina',
	       'SJAP'  => 'Taphrinamycotina',
	       'SCRY'  => 'Taphrinamycotina',
	       'SPOM'  => 'Taphrinamycotina',
	       'PJIR'  => 'Taphrinamycotina',
g	       'NIRR'  => 'Taphrinamycotina',
	       'TDEF'  => 'Taphrinamycotina',
	       'PCON'  => 'Pezizomycotina',
	       'TMEL'  => 'Pezizomycotina',
	       'ANID'  => 'Pezizomycotina',
	       'NCRA'  => 'Pezizomycotina',
	       'SCER'  => 'Saccharomycotina',
	       'CALB'  => 'Saccharomycotina',
	       'YLIP'  => 'Saccharomycotina',
	       'CCIN'  => 'Basidiomycota',
	       'PGRA'  => 'Basidiomycota',
	       'CNEO'  => 'Basidiomycota',
	       'UMAY'  => 'Basidiomycota',
	       'SROS'  => 'Basidiomycota',
	       'MBRE'  => 'Metazoa',
	       'HSAP'  => 'Metazoa');
	       
my $dir = shift || ".";

opendir(DIR,$dir) || die $!;
for my $file ( readdir(DIR) ) {
    next unless ($file =~ /(\S+)\.(gene_patterns|patterns)\.tab$/);
    my ($sp,$type) = ($1,$2);
    open(my $in => "$dir/$file" ) || die $!;
    if( $type eq 'patterns' ) {
	my %patterns;
	while(<$in>) {
	    my ($n, $pattern) = split;
	    if( $pattern ) {
		my %newpat;
		for my $t ( split(/,/,$pattern) ) {
		    if( $t ne $sp) {
			$newpat{$groups{$t} || $t}++;
		    }
		}
		$pattern = join(",", sort keys %newpat);
	    } else {
		$pattern = '';
	    }
	    $patterns{$pattern} += $n;
	}
	open(my $out => ">$dir/$sp.patterns_simple.tab" ) || die $!;
	for my $p ( sort { $patterns{$b} <=> $patterns{$a} } keys %patterns) {
	    print $out join("\t", $patterns{$p}, $p), "\n";
	}
    } elsif( $type eq 'gene_patterns' ) {
	my @genes;
	while(<$in>) {
	    my ($gene, $pattern) = split;
	    if( $pattern ) {
		my %newpat;
		for my $t ( split(/,/,$pattern) ) {
		    if( $t ne $sp ) {
			$newpat{$groups{$t} || $t}++;
		    }
		}
		$pattern = join(",", sort keys %newpat);
	    } else {
		$pattern = '';
	    }
	    push @genes, [ $gene, $pattern];
	}
	open(my $out => ">$dir/$sp.gene_patterns_simple.tab" ) || die $!;
	for my $gene ( sort { $a->[1] cmp $b->[1] } @genes ) { 
	    print $out join("\t", @$gene), "\n";
	}
    }

}
