#!/usr/bin/perl
use strict;
use warnings;

# $Id$

=head1 NAME

RIP index calculation

=head1 USAGE

RIP_index_calculation genome.fa > genome.RIP_stats.dat

=head1 DESCRIPTION

This script calculates RIP index for windows across a genome
 
RIP index types

   o "RIP Product Index" (TpA / ApT) (Margolin et al 1998, Selker et
     al 2003) 
     Sequences that have not been subjected to RIP exhibit
     values less than 0.8 while values of 1.1 indiate RIP'd sequence.

   o "RIP Substrate Index" (CpA + TpG / ApC + GpT) (Margolin et al
     1998, Selker et al 2003) 
     Sequences that have not been subject to RIP have values greater 
     than 1.1 while values less than 0.9 indicate RIP'd sequence.

   o "Composite Index" (CRI) = [product - substrate] (Lewis et al
      2009) "A positive CRI value implies that the DNA has been
      subjected to RIP" (Lewis et al).


=head1 OPTIONS

Running options:	   'r|report:s'    => \$rformat,

 -t or --type    RIP index type   ('product','substrate','composite' or CRI)
 -w or --window  Window Size       integer value, default is 500 bp
 -s or --step    Window step size  integer value < windowsize, default 100
 -m or --minlen  min RIP feat len  Minimum length of a RIP region to be reported, 
                                   can't be less than 'step' long
                                   
Input/Output
 -i or --input   Filename          Genome sequence file {REQUIRED]
 -f or --format  Sequence format   File format for BioPerl (default fasta)
 -o or --output  Output filename   Where to write the results, if 
                                   omitted will use STDOUT
 -r or --report  Report format     ('bed','wig','RIPbed','RIPgff') see OUTPUT for info

=head1 OUTPUT

There are three output formats supported.

=head2 bed

This is just a tab-delimited output in the UCSC BED format with a 
feature for each window (which are overlapping) 
http://genome.ucsc.edu/FAQ/FAQformat#format1
Score column indicates the RIP score for the reported RIP index type

=head2 wig

The UCSC WIG format http://genome.ucsc.edu/FAQ/FAQformat#format6 for
with a feature for each window (which are overlapping).

=head2 RIPbed

UCSC BED format of the identified RIP'd regions

=head2 RIPgff

GFF3 of the identified RIP'd regions


=head1 REFERENCES

Lewis ZA, Honda S, Khlafallah TK, Jeffress JK, Freitag M, Mohn F,
Schübeler D, Selker EU. Relics of repeat-induced point mutation direct
heterochromatin formation in Neurospora crassa. Genome Res 2009
19(3):427-37.

Margolin BS, Garrett-Engele PW, Stevens JN, Fritz DY, Garrett-Engele
C, Metzenberg RL, Selker EU. (1998) A methylated Neurospora 5S rRNA
pseudogene contains a transposable element inactivated by
repeat-induced point mutation. Genetics 149:1787–1797.

Selker EU, Tountas NA, Cross SH, Margolin BS, Murphy JG, Bird AP,
Freitag M. (2003) The methylated component of the Neurospora crassa
genome. Nature 422:893–897.

=head1 AUTHORS

Jason Stajich, jason.stajich [AT] ucr.edu


=cut

# BioPerl Modules
use Bio::SeqIO;
 
# System Modules
use Getopt::Long;
use List::Util qw(sum min max);

# assume zcat or bzcat are on the system of bzip of gzip files are provided
# or you may want to file these or uncompress the files
my %compress = ('gz'  => 'zcat',
		'bz2' => 'bzcat');

use constant DOUBLET_LENGTH => 2;

# variables to be set by cmdline options
my ($index_type) = qw(CRI);  # one of 'product','substrate', 'CRI' or 'composite'

my ($window,$step) = (500,10);
my $minlen         = $step;
my $sformat = 'fasta';  # sequence format
my $rformat = 'RIPgff'; # RIP report format
my %cutoffs = (composite => sub { (shift @_) >0 },
	       product   => sub { (shift @_) > 1.1},
	       substrate => sub { (shift @_) < 0.9 } );
#GFF fields
my $src;
my $type = 'region';

my (@input,$output);

# get the cmdline options
my $debug =0;
GetOptions(
	   'v|verbose|debug!' => \$debug,
	   'h|help'        => sub { exec( 'perldoc', $0);
				    exit(0); },
	   't|type:s'      => \$index_type,
	   'w|window:i'    => \$window,
	   's|step:i'      => \$step,
	   'm|minlen:i'    => \$minlen,
	   
	   'f|format:s'    => \$sformat,
	   'i|input=s'     => \@input,
	   'o|output:s'    => \$output,
	   'sourcetag:s'   => \$src,

	   'r|report:s'    => \$rformat,
	   
	   );


if( $minlen < $step) {
    warn("Minimum length cannot be less than '$step' (stepsize) long\n");
    $minlen = $step;
}
 
if( $index_type =~ /prod/i ) {
    $index_type = 'product';
} elsif( $index_type =~ /subs/i ) {
    $index_type = 'substrate';
} elsif( $index_type =~ /comp|CRI/i ) {
    $index_type = 'composite';
}

$src ||= "RIP_".$index_type."_index";

unless( @input ) {
    unless( @ARGV ) {
	die("Require an input file via -i or --input or on the cmdline. See usage with -h\n");
	
    } else {
	@input = @ARGV;
    }
}

if( $rformat !~ /wig|bed|RIPbed|RIPgff/ ) {
    die("rformat of $rformat is unrecognized. See usage with -h\n");
    
}
my $outfh;
if( defined $output ) {
    open($outfh => ">$output") || die $!;
} else {
    $outfh = \*STDOUT;
}

if($rformat eq 'RIPgff' ) {
    print $outfh "##gff-version 3\n";
}
my @doublets;
{ # some initialization code 
    for my $base1 ( qw(C A G T) ) {
	for my $base2 ( qw( C A G T) ) {
	    push @doublets, $base1.$base2;
	}
    }
}

if( $debug ) {
    warn("using $index_type indextype\n");
}

for my $file ( @input ) {
    my $fh;
    if( $debug ) {
	warn("processing $file\n");
    }
    if( $file =~ /\.(gz|bz2)$/ ) {
	open($fh, "$compress{$1} $file |") || die "Cannot open $compress{$1} $file: $!";
    } else {
	open($fh, "< $file" ) || die "Cannot open $file: $!";
    }
    if( $rformat eq 'wig' ) {
	printf $outfh "type track=%s name=\"%s\" description=\"%s %s %d bp windows\"\n",
	'wiggle_0', "$index_type RIP index",
	"$index_type",
	'RIP index calculation',
	$window;
    }
    # could change this to use Bio::DB::Fasta
    # if we require a FASTA input file
    # then wouldn't bring the whole chromosome into memory
    # not sure if it would be faster?
    my $in = Bio::SeqIO->new(-format => $sformat, 
			     -fh     => $fh);
    while( my $seq = $in->next_seq ) {	
	my @rip_features;
	my $seqid = $seq->display_id;
	warn("processing $seqid sequence\n") if $debug;
	if( $rformat eq 'wig' ) {
	    printf $outfh "fixedStep chrom=%s start=%d step=%d\n",
	    $seqid, 1, $step;
	}
	my $seqstr = $seq->seq;
	
	my $i = 0;
	my $j = 0;
	my $len = $seq->length;
	my $last_rip_start;
	while( $i <= $len ) {
	    my %freq;
	    for my $d ( @doublets ) { 
		$freq{$d} = 0;
	    }
	    my $r;
	    if( $i > ( $len - $window) ){ 
		# do some wrap around for the end-case
		#          900        1000
		# -------------------->
		#          |---------|
		#        |--------]
		#             |-------]
		#                |----]
		# FIX ME
		#
		last;
#		$r = substr($seqstr, $i-$stepsize, $window);
		
	    } else {
		$r = substr($seqstr,$i, $window);
	    }
	    my $wlen = length($r);
	    
	    for( my $k = 0; $k < $wlen - DOUBLET_LENGTH; $k++) {
		my $w = uc(substr($r,$k,DOUBLET_LENGTH));
		$freq{$w}++;
	    }
	    my %index = ( 'product'  => ( $freq{AT} ? # if not 0
					 $freq{TA} / $freq{AT} : -1),

			  'substrate'=> ( ($freq{AC} + $freq{GT} ) ? 
					  ($freq{CA} + $freq{TG}) /
					  ($freq{AC} + $freq{GT} ) : -1),
			  );
	    $index{composite} = $index{product} - $index{substrate};
	    
	    if( $rformat eq 'bed' ) {
		print $outfh join("\t",
				  $seqid,
				  $i+1,
				  $i+$window,
				  '',
				  $index{$index_type},
				  '.'),
		"\n";		
	    } elsif( $rformat eq 'wig' ) {
		printf $outfh "%.4f\n",$index{$index_type};
	    } 
	    if( $cutoffs{$index_type}->($index{$index_type}) ) {
		warn("RIP region ($index_type) has score of $index{$index_type} ($i).\n") if $debug;
		push @rip_features, [$i+1,$i+$window,
				     $index{$index_type}];
	    }
	    $i+= $step;
	}
		
	$i = 0;
	for my $f ( &collapse_rip_features(\@rip_features,$step) ) {
	    $f->[2] = sprintf("%.4f",$f->[2]); # signif digits
	    if( $rformat eq 'RIPbed' ) {
		print $outfh join("\t", 
				  $seqid,
				  $f->[0], # start
				  $f->[1], # stop
				  sprintf("%s_RIP%05d",
					  $seqid,
					  $i++),
				  $f->[2], # score
				  '.',
				  ),"\n";
	    } elsif ($rformat eq 'RIPgff' ) {
		my $ripid = sprintf("%s_RIP%05d",$seqid,$i++);
		print $outfh join("\t", $seqid,
				  $src,
				  $type,
				  @$f,
				  '.',
				  '.',
				  sprintf("ID=%s;Name=%s",
					  $ripid,
					  $ripid),
				  ), "\n";
	    }
	}	
	last if $debug;
    }
    close($fh);
    $fh = undef;
}

sub range_compare {
    # returns -1 if left < right
    # returns 0 if they overlap
    # returns 1 if left > right
    my ($left,$right) = @_;
    my @left = @$left;
    my @right = @$right;
#    $left[0] -= 1;
#    $left[1] += 1;
    
#    $right[0] -= 1;
#    $right[1] += 1;
    warn("comparing @$left to @$right\n") if $debug;
    # start of the left occurs AFTER the end of the right, can't overlap
    if( $left[0] > $right[1] ) {
	warn("returning 1 (left is after right)\n") if $debug; 
	return 1;
	# end of the left finishes BEFORE the start of the right, can't overlap
    } elsif( $left[1] < $right[0] ) {
	warn("returning -1 (left is before right)\n") if $debug; 
	return -1
    } else { 
	warn("left overlaps right\n") if $debug;
	return 0;
    }
}

sub merge { 
    my ($left,$right) = @_;
    $left->[0] = min($left->[0], $right->[0]);
    $left->[1] = max($left->[1], $right->[1]);
    push @{$left->[2]}, @{$right->[2]};
}

sub collapse_rip_features {
    my ($feats,$step_size) = @_;
    my $last;
    #               start    stop      score
    my @f = ( map { [ $_->[0], $_->[1], [ $_->[2]] ] } 
	      sort { $a->[0] <=> $b->[0] } @$feats);
    my $chg = 1;
    my $iter = 0;
    while( $chg ) {
	$chg = 0;
	for( my $i = 0; $i < scalar @f -1; $i++ ) {
	    next unless defined $f[$i];
	    for( my $j = $i+1; $j < scalar @f; $j++ ) {
		next unless defined $f[$j];
		my $cmp = &range_compare($f[$i],$f[$j]);
		if( $cmp == 0 ) {
		    &merge($f[$i],$f[$j]);
		    $f[$j] = undef;
		    $chg = 1;
		} elsif( $cmp > 0 ) {
		    last;
		}
	    }
	}
	@f = grep { defined $_ } @f;
    }

    for my $f ( @f ) {
	$f->[2] = mean(@{$f->[2]});
    }
    @f;
}

sub mean {
    return sum(@_) / scalar @_;
}

