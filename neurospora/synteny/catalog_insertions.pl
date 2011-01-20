#!/usr/bin/perl -w
use strict;
use Bio::AlignIO;
use List::Util qw(sum);
use Getopt::Long;

# pass in map file and reference genome number (order) 
# and print out indels in the reference relative to each of the other genomes.
my $min_run_length = 100;
my $debug = 0;
my $allsizes = 'allsizes.dat';
my $blockdir = 'ref_only_blocks';
GetOptions('m|min:i'    => \$min_run_length,
	   'a|all:s'    => \$allsizes,
	   'b|blockdir:s'=> \$blockdir,
	   'v|verbose!' => \$debug,);
open(my $allsizesfh => ">$allsizes") || die $!;

mkdir($blockdir);
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
my $runseen = 0;
my $blk_ct;
for my $block ( @blocks ) {    
    my $cmd = sprintf("%s %s %s %s %d %d %s",
		      'sliceAlignment', $alndir, $genomes[$ref],
		      $block->[$ref]->{chrom},
		      $block->[$ref]->{start}, $block->[$ref]->{end},
		      $block->[$ref]->{strand});
    open(my $slice => "$cmd |") || die "$cmd: $!";
    my $alnio = Bio::AlignIO->new(-fh     => $slice,
				  -alphabet => 'dna',
				  -format => 'fasta');
    if( my $aln = $alnio->next_aln ) {
	my @seqs = $aln->each_seq;
	my ($refseq);
	for my $s ( @seqs ) {
	    if( $s->id eq $genomes[$ref]) {
		$refseq = $s;
		last;
	    }
	}
	warn("refseq is ",$refseq->id, "\n") if $debug;
	my $cols = $aln->gap_col_matrix;
	my $coln = 1;
	my @ref_insert_cols;
	for my $col ( @$cols ) {
	    if( $col->{$genomes[$ref]} == 0 &&
		sum(map { $col->{$genomes[$_]} } @not_ref) == @not_ref ) {
		push @ref_insert_cols, $coln;
	    }
	    $coln++;
	}
	my @keep;
	my $offset = $block->[$ref]->{start}-1;
	for my $run ( collapse_nums(@ref_insert_cols) ) {
	    next if $run !~ /(\d+)\-(\d+)/;
	    my ($start,$end) = ($1,$2);
	    my $len = abs($end-$start);
	    if( $len >= $min_run_length ) {
		my $start_chrom_loc = $refseq->location_from_column($start);
		my $end_chrom_loc = $refseq->location_from_column($end);

		# need to convert back to CHROM coordinates		
		print $block->[$ref]->{chrom}, ": ", 
		$start_chrom_loc->end + $offset, "..", 
		$end_chrom_loc->start + $offset," length=$len\n";
		#" ($run offset=$offset - $cmd) ",
		#$start_chrom_loc->to_FTstring(), " ", $end_chrom_loc->to_FTstring(),"\n"; 
		$runseen = 1;

		my $slice_extract = sprintf("%s %s %s %s %d %d %s > %s/blk_%s.fasaln",
					    'sliceAlignment', $alndir, $genomes[$ref],
					    $block->[$ref]->{chrom},
					    $start_chrom_loc->end + $offset,
					    $end_chrom_loc->start + $offset,
					    1,$blockdir,$blk_ct++,
					    );
		`$slice_extract`;
	    }
	    print $allsizesfh $len, "\n";
	    last if $runseen && $debug;
	}
    }
    last if $runseen && $debug;
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

