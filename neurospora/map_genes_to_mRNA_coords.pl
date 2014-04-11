#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Bio::SeqIO;
use Bio::DB::Fasta;

my $gtf_file = 'neurospora_crassa_or74a_12_transcripts.gtf';
my $gff3_file = 'neurospora_crassa_or74a_12_transcripts.gff3';

my $mRNA_broad_url = 'https://www.broadinstitute.org/annotation/genome/neurospora/download/?sp=EATranscriptsGtf&sp=SNC12&sp=S.zip';

my $genome_broad_url = 'http://www.broadinstitute.org/annotation/genome/neurospora/download/?sp=EASupercontigsFasta&sp=SNC12&sp=S.zip';
my $transcriptfasta_broad_url = 'http://www.broadinstitute.org/annotation/genome/neurospora/download/?sp=EATranscriptsFasta&sp=SNC12&sp=S.zip';

my $outfile = 'exons_in_mRNA_coords.bed';
my $outfasta = 'mRNA.fasta';
my $genome = 'neurospora_crassa_or74a_12_supercontigs.fasta';

my $debug = 0;
my $skip_alt = 1;
GetOptions('g|gtf:s'     => \$gtf_file,
	   'gff3:s'      => \$gff3_file,
	   'url:s'       => \$mRNA_broad_url,
	   'o|out:s'     => \$outfile,
	   'db|genome:s' => \$genome,
	   'f|fasta:s'   => \$outfasta,
	   'v|debug!'    => \$debug,
	   'skipalt!'    => \$skip_alt,
	   );

if ( ! -f $gtf_file ) {
  `curl "$mRNA_broad_url" > gtf.zip`;
  `unzip gtf.zip`;
}
if ( ! -f $gff3_file ) {
  `perl gtf2gff3_3level.pl $gtf_file > $gff3_file`;
}
my $dbh = Bio::DB::Fasta->new($genome);

my $outseq = Bio::SeqIO->new(-format => 'fasta',
			     -file => ">$outfasta");

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
    $genes{$id}->{seqid} = $row[0];
  } elsif ( $row[2] eq 'mRNA') {
    $transcripts{$id}->{parent} = $parent;
    $transcripts{$id}->{seqid} = $row[0];
    $transcripts{$id}->{start} = $row[3];
    $transcripts{$id}->{end} = $row[4];
    $transcripts{$id}->{strand} = $row[6] eq '+' ? 1 : -1;
  } else {
    push @{$transcripts{$parent}->{'children'}->{$row[2]}}, { type  => $row[2],
							      seqid => $row[0],
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
  next if $skip_alt && $transcriptid !~ /T0$/;
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
  my $transcript_seq;
  foreach my $f ( sort { $a->{start} * $a->{strand} <=> $b->{start} * $b->{strand} } @all_features ) {
    my $start = $f->{start};
    my $end   = $f->{end};
    my $exonlen = $end - $start;
    warn($start, "..",$end, " ", $f->{strand}, " ", $transcriptid,"\n");
    if ( $f->{strand} < 0 ) {
      ($start,$end) = ($end,$start);
    }
    my $exon = $dbh->seq($f->{seqid}, $start => $end);
    $transcript_seq .= $exon;
    print $ofh join("\t", # $start,$end,
	       $transcriptid,
	       $running, ($running+$exonlen),sprintf("%s_%s_%d",
						     $transcriptid,
						     $f->{type},
						     ++$typect{$f->{type}}),
		   ),"\n";
    $running += $exonlen+1;
  }
  $outseq->write_seq(Bio::PrimarySeq->new(-id => $transcriptid,
					  -seq => $transcript_seq));
  last if $debug;
}

__END__
