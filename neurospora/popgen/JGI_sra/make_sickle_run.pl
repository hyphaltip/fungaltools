#!/usr/bin/perl -w
use strict;
use File::Spec;
my $odir = File::Spec->rel2abs('../fastq');
my $dir = shift || ".";
opendir(DIR,$dir) || die $!;
for my $d2 ( readdir(DIR) ) {
    next unless $d2 =~ /^FGSC/;
    opendir(D2, "$dir/$d2") || die $!;
    my %files;
    for my $file ( readdir(D2) ) {
	next unless ($file =~ /(\S+)\.fastq\.gz/);
	my $stem = $1;
	if( $stem =~ s/(\S+)\_(\d+)// ) {
	    $files{$1}->{pe}->{$2} = $file;
	} else {
	    $files{$stem}->{se} = $file;
	}
    }
    my $fulldir = File::Spec->rel2abs("$dir/$d2");
    open(my $ofh => ">$dir/$d2/run_sickle.sh") || die $!;
    print $ofh "#PBS -l nodes=1:ppn=1 -d $fulldir -q js -N $d2.sickle\n";
    print $ofh "module load sickle\n";

    while ( my ($stem,$d) = each %files ) {
	if( exists $d->{se} ) {
	    printf $ofh "sickle se -f %s -o $odir/%s.unpaired_se.fq -t sanger -l 25 -q 30\n",$d->{se}, $stem;
	} 
	if( exists $d->{pe} ) {
	    printf $ofh "sickle pe -f %s -r %s -o $odir/%s_%s_1.fq -p $odir/%s_%s_2.fq -s $odir/%s.unpaired_pe.fq -t sanger -l 25 -q 30\n",
	    $d->{pe}->{1}, $d->{pe}->{2},
	    $d2,$stem, $d2,$stem,$stem;
	}
	printf $ofh "cat $odir/%s.unpaired_se.fq $odir/%s.unpaired_pe.fq > $odir/%s_%s.single.fq\n",$stem,$stem,$d2,$stem;
	printf $ofh "rm $odir/%s.unpaired_se.fq $odir/%s.unpaired_pe.fq\n",
	$stem,$stem;
    }
    
}
