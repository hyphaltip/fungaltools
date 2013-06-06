#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
my $agp = shift;
my $fasta = shift;

my $in = Bio::SeqIO->new(-format => 'fasta', -file => $fasta);
my $out = Bio::SeqIO->new(-format => 'fasta');

open(my $fh => $agp) || die $!;
my %keep;
while(<$fh>) {
    next if /^\#/;
    my @row = split;
    next if $row[4] eq 'N';
    $keep{$row[5]}++;
}

while( my $seq = $in->next_seq ) {
    next unless $keep{$seq->id};
    $out->write_seq($seq);
}
