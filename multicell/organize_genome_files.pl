#!/usr/bin/perl -w
use strict;
use Bio::DB::Taxonomy;
use Getopt::Long;
use Bio::SeqIO;
use File::Spec;

my %known = ('Piromyces sp' => 73868,
	     'Nematocida sp' => 944018,
	     'Salpingoeca sp' => 946362);

my $taxonfiledir = '/shared/stajichlab/db/taxonomy';
my $indexdir;
my $final_dir;
GetOptions(
    'taxondir:s'  => \$taxonfiledir,
    'index:s'     => \$indexdir,
    'f|target|final:s' => \$final_dir,
    );

my $dir = shift || die $!;
$final_dir = shift @ARGV unless defined $final_dir;
die "need a target dir" unless $final_dir;
#my $remotedb = Bio::DB::Taxonomy->new(-source => 'entrez');
my $taxdb = Bio::DB::Taxonomy->new(-source => 'flatfile',
                                   -nodesfile => File::Spec->catfile
                                   ($taxonfiledir,'nodes.dmp'),
                                   -namesfile => File::Spec->catfile
                                   ($taxonfiledir,'names.dmp'),
                                   -directory => $indexdir || $taxonfiledir,
                                   );


mkdir($final_dir) unless -d $final_dir;

opendir(DIR => $dir) || die $!;

print join("\t", qw(SUPERKINGDOM KINGDOM PHYLUM SUBPHYLUM GENUS SPECIES PREFIX TAXONID PEPCOUNT)),"\n";
for my $file ( readdir(DIR) ) {
    next unless $file =~ /(\S+)\.fa(s|sta)?$/;
    my $stem = $1;
    my ($genus,$species,$rest) = split(/[\._]/,$stem,3);

    my $ncbi_taxid = $taxdb->get_taxonid("$genus $species");
    if( ! $ncbi_taxid && $known{"$genus $species"} ) {
	 $ncbi_taxid = $known{"$genus $species"};
    }
    my $taxnode = $taxdb->get_Taxonomy_Node(-taxonid => $ncbi_taxid);
    if( ! $taxnode ) {
	$ncbi_taxid = $taxdb->get_taxonid("$genus $species $rest");
	if( ! $ncbi_taxid) {
	 warn("Cannot find taxon node for $genus $species ($file)\n");
	 next;
	} else {
	 $taxnode = $taxdb->get_Taxonomy_Node(-taxonid => $ncbi_taxid);
	}
    } 
    my $tree = Bio::Tree::Tree->new(-node => $taxnode);
    my @lineage = $tree->get_lineage_nodes($taxnode);
    my ($speciesnode) = @{$tree->find_node(-rank => 'species')->name('scientific')};
	$speciesnode =~ s/[\"']//g;
    my @node_names = split(/\s+/,$speciesnode,2);
    my $prefix = sprintf("%s%s",uc substr($node_names[0],0,1),
			 lc substr($node_names[1],0,3));
    shift @lineage;
    shift @lineage;
    unshift @node_names, (map { my ($val) = @{$_->name('scientific')};
			   $val; } $lineage[0], $lineage[1], $lineage[2]);
    if( $node_names[2] eq 'Dikarya' ) {
       splice(@node_names, 3, 0,$lineage[3]->scientific_name);
    } else {
      splice(@node_names, 3, 0,'');
    }

#( ( map { my $node = $tree->find_node(-rank => $_);
#				my $val = 'NONE';
#				if( $node ) {
#				    ($val) = @{$node->name('scientific')};
#				}
#				$val } 
#			 qw(kingdom phylum)));
    
    my $count;
    my $out = Bio::SeqIO->new(-format => 'fasta',
			     -file   => sprintf(">%s/%s_%s.proteins.fasta",
						$final_dir,
						$node_names[-2],
						$node_names[-1]));
    my $in = Bio::SeqIO->new(-format => 'fasta',
			     -file   => File::Spec->catfile($dir,$file));
    while(my $seq = $in->next_seq) {
	my $id = $seq->display_id;
	$id =~ s/:pep//;
	if( $id =~ /^([^\|]+)\|(\S+)/) {
	    $id = $2;
	}

	$seq->display_id(sprintf("%s|%s",$prefix,$id));
	$out->write_seq($seq);
	$count++;
	#last;
    }
    print join("\t",@node_names,$prefix,$ncbi_taxid,$count),"\n";
}



