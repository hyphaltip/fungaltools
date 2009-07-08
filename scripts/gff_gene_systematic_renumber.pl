#!/usr/bin/perl -w

=head1 NAME

gff_gene_systematic_renumber

=head1 USAGE

 -l or -t or --locus or --tag    Locus Tag prefix (i.e. NCU)
 -d or --db                      Fasta genome database for getting contig lens
 -o or --output                  [Optional] output filename (STDOUT is default)
 --offset                        Starting number of the gene numbering counter

=head1 DESCRIPTION

=head1 AUTHOR

Jason Stajich jason.stajich [at] ucr.edu

=cut

use strict;
# system modules
use Getopt::Long;
# bioperl modules
use Bio::DB::Fasta;

my $locus_tag;
my $fasta_db;
my $counter_offset = 0;
my $output;
GetOptions(
	   'l|t|locus|tag:s' => \$locus_tag,
	   'd|db:s'          => \$fasta_db,
	   'o|output:s'      => \$output,
	   'offset:i'        => \$counter_offset,
	  ); 

die 'must provide a locus tag prefix with -l or --locus' unless defined $locus_tag;
die 'must provide a fastadb with -d or --db' unless defined $fasta_db;

my $db = Bio::DB::Fasta->new($fasta_db);

my ($genename,%order,%genes);
while(<>) {
    next if /^\#/;
    chomp;
    my @line = split(/\t/,$_);
    my %lastcol = map { split(/=/,$_) } split(/;/,pop @line);

    my $id = $lastcol{'ID'};    
    die("No ID for line $_\n") unless defined $id;

    if( $line[2] eq 'gene' ) {
	$genename = $id;
	push @{$order{$line[0]}}, [$line[3],$genename];
    }
    my $name =  $lastcol{'Name'} || $lastcol{'Parent'};
    my @rest = sort grep { ! /^ID|Name|Parent$/ } keys %lastcol;   
    my $lastcol = join(";", map { sprintf("%s=%s",$_,$lastcol{$_}) } @rest);

    push @{$genes{$genename}->{$line[2]}}, [@line,$id,$name,$lastcol];
    
}

my @chroms = sort { $db->length($b) <=> $db->length($a) } keys %order;

my $ofh;
if( $output ) {
    open($ofh, ">$output") || die $!;
} else {
    $ofh = \*STDOUT;
}
my $gcounter = $counter_offset;
my $tcounter = $counter_offset;
print $ofh "##gff-version 3\n";
my %old2new_geneid;
my %old2new_tid;
for my $chrom ( @chroms ) {
    for my $gene ( sort { $a->[0] <=> $b->[0] } @{$order{$chrom}} ) {
	my ($start,$genename) = @$gene;
	my ($generow) = @{$genes{$genename}->{gene}};
	my $gextra   = pop @$generow;
	my $gname    = pop @$generow;
	my $gid      = pop @$generow;
	$gname = sprintf("%sG_%05d",$locus_tag, $gcounter++);
	$old2new_geneid{$gid} = $gname;
	my $lastcol = sprintf("ID=%s;Name=%s",$gname,$gname);
	$lastcol .= ";$gextra" if defined $gextra && length $gextra;
	print $ofh join("\t", @$generow,$lastcol),"\n";
	for my $mRNA ( @{$genes{$genename}->{'mRNA'}} ) {
	    my $textra   = pop @$mRNA;
	    my $tname    = pop @$mRNA;
	    my $tid      = pop @$mRNA;
	    $tname = sprintf("%sT_%05d",$locus_tag, $tcounter++);
	    $old2new_tid{$tid} = $tname;
	    my $lastcol = sprintf("ID=%s;Parent=%s;Name=%s",
				  $tname,$gname,$tname);
	    $lastcol .= ";$textra" if defined $textra && length $textra;
	    print $ofh join("\t", @$mRNA,$lastcol),"\n";
	}
	for my $e ( sort {$a->[3] <=> $b->[3] } 
		       @{$genes{$genename}->{'CDS'}}, 
		       @{$genes{$genename}->{'exon'}} ) {
	    my $eextra   = pop @$e;
	    my $ename    = pop @$e;
	    my $eid      = pop @$e;
	    my $lastcol = sprintf("Parent=%s",
				  $old2new_tid{$ename});
	    $lastcol .= ";$eextra" if defined $eextra && length $eextra;
	    print $ofh join("\t", @$e,$lastcol),"\n";
	}
    }
}

