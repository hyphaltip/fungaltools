#!/usr/bin/perl -w
use strict;

my $datfile = shift || die $!;
my $outbase = $datfile;
$outbase =~ s/\.(tsv|tab)$//;

my @timepoints = (0,20,40,60,80,100,120,140);

open(my $fh => $datfile) || die $!;
my $header =<$fh>;

my %groups;
while (<$fh>) {
  my ($gene,@data) = split;
  my @cols = find_max_col(@data);
  my $peak1 = $cols[0]->[0];
  my $peak2 = $cols[1]->[1];

  push @{$groups{$timepoints[$peak1]}}, [$gene,@data];
}

for my $group ( sort { $a <=> $b } keys %groups ) {
  open(my $ofh => ">$outbase.peak_$group.tab") || die $!;
  print $ofh join("\t", qw(GENE), map { sprintf("TIME_%d",$_) } @timepoints),
    "\n";
  for my $d ( @{$groups{$group}} ) {
    print $ofh join("\t", @$d), "\n";
  }
}

sub find_max_col {
  my $n = 0;
  # data structure is sorted by value and each val
  # is a pair of the column # for the max and the value
  my @dat = sort { $b->[1] <=> $a->[1] } map { [ $n++, $_] } @_;
  return @dat;
}
