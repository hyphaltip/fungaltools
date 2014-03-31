#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my $gtf_file = 'neurospora_crassa_or74a_12_transcripts.gtf';
my $gff3_file = 'neurospora_crassa_or74a_12_transcripts.gff3';

my $mRNA_broad_url = 'https://www.broadinstitute.org/annotation/genome/neurospora/download/?sp=EATranscriptsGtf&sp=SNC12&sp=S.zip';

my $outfile = 'exons_in_mRNA_coords.bed';

GetOptions('g|gtf:s' => \$gtf_file,
	   'gff3:s'  => \$gff3_file,
	   'url:s'   => \$mRNA_broad_url,
	   'o|out:s' => \$outfile,
	   );

if ( ! -f $gtf_file ) {
  `curl "$mRNA_broad_url" > gtf.zip`;
  `unzip gtf.zip`;
}
if ( ! -f $gff3_file ) {
  `perl gtf2gff3_3level.pl $gtf_file > $gff3_file`;
}

open(my $ofh => ">$outfile") || die "cannot open $outfile: $!";


open(my $fh => $gff3_file ) || die "cannot open $gff3_file: $!";

my (%genes,%transcripts);
my $i = 1;
while (<$fh>) {
  next if /^\#/;
  chomp;
  my @row = split(/\t/,$_);
  my $last_col = pop @row;
  my %groups = map { split(/=/,$_) } split(/;/,$last_col);
  my ($id,$parent) = ($groups{'ID'},$groups{'Parent'});
  if ( ! defined $id && ! defined $parent ) {
    warn("feature doesn't have a id or parent line [$i] $_\n");
    next;
  }
  if ( $row[2] eq 'gene' ) {
    my $name  = $groups{'Name'};
    if ( $name ne $id ) {
      warn("$name ne $id\n");
    }
    $genes{$id}->{start} = $row[3];
    $genes{$id}->{end} = $row[4];
    $genes{$id}->{strand} = $row[6] eq '+' ? 1 : -1;
    $genes{$id}->{name} = $name;
  } elsif ( $row[2] eq 'mRNA') {
    $transcripts{$id}->{parent} = $parent;
    $transcripts{$id}->{start} = $row[3];
    $transcripts{$id}->{end} = $row[4];
    $transcripts{$id}->{strand} = $row[6] eq '+' ? 1 : -1;
  } else {
    push @{$transcripts{$parent}->{'children'}->{$row[2]}}, { type  => $row[2],
							       start => $row[3],
							       end   => $row[4],
							       strand=> $row[6] eq '+' ? 1 : -1,
							       phase => $row[7]};
  }
  $i++;
}

warn "processed $i lines\n";

while ( my ($transcriptid,$hash) = each %transcripts ) {
#  print "$transcriptid\n";
  my @all_features;
  
  while ( my ($ftype,$features) = each %{$hash->{children}} ) {
#    printf "%s has %d features for %s\n",$ftype,scalar @{$features},
#	$transcriptid;
    next if $ftype eq 'exon';
    push @all_features, @{$features};
  }
  my $offset = $hash->{start};
  my $running = 1;
  my %typect;
  foreach my $f ( sort { $a->{start} * $a->{strand} <=> $b->{start} * $b->{strand} } @all_features ) {
    my $start = $f->{start};
    my $end   = $f->{end};
    my $exonlen = $end - $start;
    print $ofh join("\t", # $start,$end, 
	       $transcriptid,
	       $running, ($running+$exonlen),sprintf("%s_%s_%d",
						     $transcriptid,
						     $f->{type},
						     ++$typect{$f->{type}}),
		   ),"\n";
    $running += $exonlen+1;
  }
}

__END__
#IGNORE
open(my $fh => $gtf_file ) || die "cannot open $gtf_file: $!";

my %genes;
my $i = 1;
while (<$fh>) {
  next if /^\#/;
  chomp;
  my $msg;
  if ( s/#\s*(.+)// ) {
    $msg = $1;
    warn("($msg) for $_\n");
  }
  my @row = split(/\t/,$_);
  my $last_col = pop @row;
  my %groups = map { split(/\s+/,$_) } split(/;\s*/,$last_col);
  if ( keys %groups > 2) {
    warn("last_col was $last_col\n");
  }
  my ($gene_id,$transcript_id) = ($groups{'gene_id'},
				  $groups{'transcript_id'});
  if ( ! defined $gene_id || ! defined $transcript_id ) {
    warn("feature doesn't have a gene or transcript id, may need to skip it, line [$i]\n");
    next;
  }
  $transcript_id =~ s/\"//g;
  $gene_id =~ s/\"//g;
  push @{$genes{$gene_id}->{$transcript_id}->{$row[2]}}, [@row];
  $i++;
}

warn "processed $i lines\n";

while ( my ($gene,$trans) = each %genes ) {
  print "$gene\n";
  for my $transid  ( keys %$trans ) {
    print "\t$transid\n";
    for my $subfeat ( sort keys %{$trans->{$transid}} ) {
      print "\t\t$subfeat\n";
    }
  }
}
