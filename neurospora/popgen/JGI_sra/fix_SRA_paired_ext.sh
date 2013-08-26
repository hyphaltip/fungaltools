perl -i -p -e 's/^([@\+]SRR\d+\.\d+)/$1\/1/' FGSC*/*_1.fastq
perl -i -p -e 's/^([@\+]SRR\d+\.\d+)/$1\/2/' FGSC*/*_2.fastq
