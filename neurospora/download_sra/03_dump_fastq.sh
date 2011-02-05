cd SRP004848
rmdir SRR*
for file in *.lite.sra
do
 base=`basename $file .lite.sra`
 fastq-dump -alt 1 -A $base $file
done
