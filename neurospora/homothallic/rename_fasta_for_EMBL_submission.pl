#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
# rename the scaffold names to be simple and consistent with what is submitted to GENBANK/EMBL
my $dir = 'mercator';
opendir(D,$dir) || die $!;
for my $f ( readdir(D) ) {
    if( $f =~ /(\S+)\.fa/) {
	my $stem = $1;
	my $in = Bio::SeqIO->new(-format => 'fasta',
				 -file   => "$dir/$f");
	my $out = Bio::SeqIO->new(-format => "fasta",
				  -file   => ">$stem.fasta");
	while(my $seq = $in->next_seq ) {
	    my $id = $seq->display_id;
	    
	    if( $id =~ /^assembled(\d+)/) {
		$id = sprintf("SC_%d",$1);
	    } elsif( $id =~ /^NODE_(\d+)/) {
		$id = sprintf("scf_%d",$1);
	    } else {
		warn("unmatched scaffold name $id\n");
	    }	    
	    $seq->display_id($id);
	    $out->write_seq($seq);
	}
    }
}
