#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Bio::DB::Fasta;

use Getopt::Long;

my ($dir,$input,$db,$debug);
GetOptions('d|dir:s' => \$dir,
	   'i|input:s' => \$input,
	   'db:s'     => \$db,
	   'v|debug!' => \$debug,
	  );


unless  ( $db ) {
  die("need a db to proceed - use -db\n");
}
my $dbh = Bio::DB::Fasta->new($db);
$input ||= shift @ARGV || die "need an input";

if ( ! $dir ) {
  ($dir) = split(/\./,$input);
}
mkdir($dir) unless -d $dir;


open(my $in => $input) || die $!;

while (<$in>) {
  my ($gene,@genes) = split;
  my %nms = map { split(/=/,$_) } @genes;
  my $pattern = join('-',map { (split(/\|/,$_))[0] } sort keys %nms);
  if ( ! -d "$dir/$pattern") {
    mkdir("$dir/$pattern");
  }
  my ($sp,$genename) = split(/\|/,$gene);
  my $out = Bio::SeqIO->new(-format => 'fasta',
			    -file   => ">$dir/$pattern/$genename.genes.fa");
  if ( ! exists $nms{$gene} ) {
    warn("did not have self ($gene) in the list for the gene?\n");
    warn("\t",join("\n\t",keys %nms),"\n");
    $nms{$gene} = 0;
  } 
  for my $genename ( keys %nms ) {
    next if $genename eq $gene;
    my $seq = $dbh->get_Seq_by_id($genename);
    if ( ! defined $seq) {
      warn("cannot find $genename\n");
    } else {
      $out->write_seq($seq);
    }
  }
  last if $debug;
}
