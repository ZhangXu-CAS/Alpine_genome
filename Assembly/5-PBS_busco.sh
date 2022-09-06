#PBS -N Sau_busco
#PBS -l nodes=1:ppn=20
#PBS -q batch 

PBS_O_WORKDIR=/public/home/WHC_zhangx/Saussurea_genome
cd $PBS_O_WORKDIR/BUSCO
CONG=$PBS_O_WORKDIR/assembly_files/NextDenovo/Saussurea_assemble2.1/03.ctg_graph
DB=$PBS_O_WORKDIR/BUSCO/database/embryophyta_odb10
time busco -m geno -i $CONG  -l $DB -o embryophyta_odb10_busco -c 20  --offline
echo "busco done"
