#PBS -N nd_Pilon
#PBS -l nodes=1:ppn=2
#PBS -q fat 

PBS_O_WORKDIR=/public/home/WHC_zhangx/Saussurea_genome/polish/Pilon
cd $PBS_O_WORKDIR
R1=/public/home/WHC_zhangx/Saussurea_genome/illumina_reads/sr.1.fq.gz
R2=/public/home/WHC_zhangx/Saussurea_genome/illumina_reads/sr.2.fq.gz

bwa index -p index/necat_racon3 index/necat_racon3.fasta
bwa mem -t 40  index/necat_racon3  $R1 $R2 |samtools view -Sb > necat_racon3.bam
samtools sort -@ 40 -O bam -o necat_racon3_sort.bam necat_racon3.bam
samtools index -@ 40 necat_racon3_sort.bam

java -Xmx500G -jar pilon-1.23.jar --genome index/necat_racon3.fasta --frags necat_racon3_sort.bam --diploid --threads 2 --output necat_racon3_pilon 
echo "Pilon done"


