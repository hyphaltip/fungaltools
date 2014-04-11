module load bedtools

bedtools coverage -a ../repeats/neurospora_crassa_or74a_12_supercontigs.fasta.RM_filter.gff3 -b raw/Nc.windows.10kb -counts > raw/Nc.repeats.10kb.counts.txt
bedtools coverage -a ../repeats/neurospora_crassa_or74a_12_supercontigs.Gioti.RM_filter.gff3 -b raw/Nc.windows.10kb -counts > raw/Nc.repeats.10kb.Gioti.counts.txt
perl /rhome/jstajich/src/genome-scripts/viz/normalize_window_readcount.pl raw/Nc.repeats.10kb.*.counts.txt
perl -i -p -e 's/^Super/Nc_Super/' raw/Nc.repeats.10kb.*.norm.txt
cp raw/Nc.repeats.10kb.*.norm.txt data/repeats
