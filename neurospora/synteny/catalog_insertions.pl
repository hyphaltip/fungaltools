#!/usr/bin/perl -w
use strict;
use Bio::AlignIO;
use List::Util qw(sum);
use Getopt::Long;

# pass in map file and reference genome number (order) 
# and print out indels in the reference relative to each of the other genomes.
my $min_run_length = 100;
open(my $allsizesfh => ">allsizes.dat") || die $!;

GetOptions('m|min:i' => \$min_run_length);

my ($alndir,$ref) = @ARGV;

$ref ||= 0;
my @genomes;

{
    open(my $g => "$alndir/genomes") || die $!;
    while(<$g>) {
	push @genomes, split;
    }
    close($g);
}

open(my $fh => "$alndir/map") || die $!;
my @blocks;
while(<$fh>){
    my ($blockid, @syntenic_blocks) = split;
    my @blockset;
    while( @syntenic_blocks ) {
	my ($chrom,$start, $end, $strand) = splice(@syntenic_blocks,0,4);
	push @blockset, { 'chrom' => $chrom,
			  'start' => $start,
			  'end'   => $end,
			  'strand'=> $strand};
    }
     push @blocks, [@blockset];
}
my @not_ref = grep { $_ != $ref } 0..$#genomes;

for my $block ( @blocks ) {
    my $cmd = sprintf("%s %s %s %s %d %d %s",
		      'sliceAlignment', $alndir, $genomes[$ref],
		      $block->[$ref]->{chrom},
		      $block->[$ref]->{start}, $block->[$ref]->{end},$block->[$ref]->{strand});
    open(my $slice => "$cmd |") || die "$cmd: $!";
    my $alnio = Bio::AlignIO->new(-fh     => $slice,
				  -alphabet => 'dna',
				  -format => 'fasta');
    if( my $aln = $alnio->next_aln ) {
	my @seqs = $aln->each_seq;
	my $cols = $aln->gap_col_matrix;
	my $coln = 0;
	my @ref_insert_cols;
	for my $col ( @$cols ) {
	    if( $col->{$genomes[$ref]} &&
		sum(map { $col->{$genomes[$_]} } @not_ref) == 0) {
		push @ref_insert_cols, $coln;
	    }
	    $coln++;
	}
	my @keep;
	for my $run ( collapse_nums(@ref_insert_cols) ) {
	    next if $run !~ /(\d+)\-(\d+)/;
	    my ($start,$end) = ($1,$2);
	    my $len = abs($end-$start);
	    if( $len >= $min_run_length ) {
		print $block->[$ref]->{chrom}, ": ", $start + $block->[$ref]->{start}, "..", $end + $block->[$ref]->{start}, "\n";
	    }
	    print $allsizesfh $len, "\n";
	}
    }
}   

# from Bio::Search::SearchUtils
sub collapse_nums {
# This is probably not the slickest connectivity algorithm, but will do for now.
    my @a = @_;
    my ($from, $to, $i, @ca, $consec);
    
    $consec = 0;
    for($i=0; $i < @a; $i++) {
        not $from and do{ $from = $a[$i]; next; };
    # pass repeated positions (gap inserts)
    next if $a[$i] == $a[$i-1];
        if($a[$i] == $a[$i-1]+1) {
            $to = $a[$i];
            $consec++;
        } else {
            if($consec == 1) { $from .= ",$to"; }
            else { $from .= $consec>1 ? "\-$to" : ""; }
            push @ca, split(',', $from);
            $from =  $a[$i];
            $consec = 0;
            $to = undef;
        }
    }
    if(defined $to) {
        if($consec == 1) { $from .= ",$to"; }
        else { $from .= $consec>1 ? "\-$to" : ""; }
    }
    push @ca, split(',', $from) if $from;

    @ca;
}

