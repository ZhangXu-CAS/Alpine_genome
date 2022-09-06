#PBS -N PCOC_DOCKER
#PBS -l nodes=1:ppn=20
#PBS -q comput

WD=/public/home/WHC_zhangx/Saussurea_genome/Alpine/pcoc
cd $WD
module load apps/singularity/gnu/2.4.2
#singularity pull docker://carinerey/pcoc
CMD_PCOC_DOCKER="singularity exec ./pcoc.simg"
#$CMD_PCOC_DOCKER  pcoc_num_tree.py -t  mydata/OG0000000.tree -o mydata/OG0000000_num_tree.pdf
scenario="1/3/8/13/22/27/33"
cat pep3_orth_names | while read OG
do $CMD_PCOC_DOCKER pcoc_det.py -t $WD/mydata/trees/$OG.tre -aa $WD/mydata/seqs/$OG.mafft.fa -o $WD/mydata/outputs/${OG}_pcoc -m $scenario -f 0.9  -cpu 40
done
echo "PCOC_DOCKER done"


#cat ../pep3_orth_names | while read name
#do cp /public/home/WHC_zhangx/Saussurea_genome/Alpine/pcoc/mydata/outputs/${name}_pcoc/RUN_*/*.filtered_results.tsv .
#done
# $CMD_PCOC_DOCKER pcoc_det.py -t mydata/trees/OG0000000.tre -aa mydata/seqs/OG0000000.mafft.fa -o mydata/OG0000000_pcoc -m "1/3/8/13/22/27/33" -f 0.9  --plot_complete_ali --plot