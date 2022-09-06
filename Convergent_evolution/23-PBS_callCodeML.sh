#PBS -N callCodeML
#PBS -l nodes=1:ppn=2
#PBS -q comput


PBS_O_WORKDIR=/public/home/WHC_zhangx/Saussurea_genome/Alpine/PAML

cd $PBS_O_WORKDIR
python3 callCodeml2.py msa_cds_pep3_2 all_alpine.tre
echo "Sau_callCodeML done"
