#!/usr/bin/perl -w
use strict;

use File::Spec;
my $segexe = `which seg`;
if( ! -x $segexe ) {
    $segexe = '/opt/wu-blast/2009-07/filter/seg';
}
my $dir = 'clean';
my $odir = 'seg';
my $ext = 'seg';
opendir(DIR,$dir)|| die $!;
for my $file ( readdir(DIR) ) {
    next unless $file = /(\S+)\.fasta$/;
    my $base = $1;
    my $ofile = File::Spec->catfile($odir,"$base.$ext");
    next if -f $ofile;
    warn "$segexe $dir/$file -q -z 1  > $ofile";
}
