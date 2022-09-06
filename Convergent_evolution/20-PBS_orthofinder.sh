#PBS -N orthofinder
#PBS -l nodes=1:ppn=40
#PBS -q comput 


PBS_O_WORKDIR=/public/home/WHC_zhangx/Saussurea_genome/Alpine

cd $PBS_O_WORKDIR
orthofinder -t 80 -a 80 -M msa -f peps3

echo "orthofinder done"
