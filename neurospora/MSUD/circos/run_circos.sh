#PBS -l nodes=1:ppn=1,mem=2gb -N circos -j oe
module load circos 
circos -conf etc/circos.conf
