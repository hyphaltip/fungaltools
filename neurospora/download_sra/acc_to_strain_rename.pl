#!/usr/bin/perl -w
use strict;
my $names = 'SRAaccession_to_strains.txt';

# CMDLINE input
my $dir = shift || "."; # expect the directory of all the SRR files to be provided

open(my $fh => $names) || die $!;
my %names;
while(<$fh>) {
 my ($accn,$strain,$other) = split(/\t/,$_);
 $strain =~ s/\s+/_/g;
 $names{$accn} = $strain;
}

opendir(DIR, $dir) || die $!;
for my $d ( readdir(DIR) ) {
 next unless -d "$dir/$d" && $names{$d};
 print "mv $dir/$d/$d.lite.sra $dir/$names{$d}.lite.sra\n";
 `mv $dir/$d/$d.lite.sra $dir/$names{$d}.lite.sra`;
}

