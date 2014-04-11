module load bedtools

#bedtools coverage -abam ../allreads_aln_Nc12/PP.bam -b raw/Nc.windows.100kb -counts > raw/Nc.PP.100kb.counts.txt
#bedtools coverage -abam ../allreads_aln_Nc12/2PF.bam -b raw/Nc.windows.100kb -counts > raw/Nc.2PF.100kb.counts.txt
#bedtools coverage -abam ../allreads_aln_Nc12/4PF.bam -b raw/Nc.windows.100kb -counts > raw/Nc.4PF.100kb.counts.txt

bedtools coverage -abam ../allreads_aln_Nc12/PP.bam -b raw/Nc.windows.10kb -counts > raw/Nc.PP.10kb.counts.txt
bedtools coverage -abam ../allreads_aln_Nc12/2PF.bam -b raw/Nc.windows.10kb -counts > raw/Nc.2PF.10kb.counts.txt
bedtools coverage -abam ../allreads_aln_Nc12/4PF.bam -b raw/Nc.windows.10kb -counts > raw/Nc.4PF.10kb.counts.txt

perl /rhome/jstajich/src/genome-scripts/viz/normalize_window_readcount.pl raw/*.counts.txt
