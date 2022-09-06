#PBS -N ipr_Rheum_alexandrae
#PBS -l nodes=1:ppn=20
#PBS -q fat 
source ~/.bashrc
PBS_O_WORKDIR=/public/home/WHC_zhangx/Saussurea_genome/annotation
cd $PBS_O_WORKDIR
./interproscan-5.51-85.0/interproscan.sh  -i Rheum_alexandrae.fa -cpu 80 -b Rheum_alexandrae.ipr -goterms -iprlookup -pa -dp
echo 'ipr done'
