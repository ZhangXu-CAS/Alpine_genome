#PBS -N Sau_GCE
#PBS -l nodes=1:ppn=40
#PBS -q comput 

PBS_O_WORKDIR=/public/home/WHC_zhangx/Saussurea_survey
cd $PBS_O_WORKDIR
./kmerfreq-master/kmerfreq  \
  -k 19 -t 10 -p ss reads_files.lib &> 19kmer_freq.log
./gce-1.0.0/gce \
  -f ss.freq.stat -c 25 -g 44080884843 -m 1 -D 8 -b 1 > ss.table 2> ss.log
 echo "GCE done"
 