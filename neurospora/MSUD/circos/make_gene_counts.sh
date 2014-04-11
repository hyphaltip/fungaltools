module load bedtools

grep gene ../Nc12/neurospora_crassa_or74a_12_transcripts.gff3 > raw/Nc12.genes.gff3
bedtools coverage -a raw/Nc12.genes.gff3 -b raw/Nc.windows.100kb -counts > raw/Nc.genes.100kb.counts.txt
#perl /rhome/jstajich/src/genome-scripts/viz/normalize_window_readcount.pl raw/Nc.genes.100kb.counts.txt
perl -i -p -e 's/^Super/Nc_Super/' raw/Nc.genes.100kb.counts.txt
cp raw/Nc.genes.100kb.counts.txt data/genes
