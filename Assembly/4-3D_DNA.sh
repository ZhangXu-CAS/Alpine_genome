#!~/bin/bash

d=/home/wangh/zx/Saussurea_genome/hic/3ddna/3d-dna
juicer=/home/wangh/zx/Saussurea_genome/hic/3ddna/juicer
draft='ndasm2_ndph_pilon_purged'
python $juicer/misc/generate_site_positions.py DpnII  ${draft}  ${draft}.fa
awk 'BEGIN{OFS="\t"}{print $1, $NF}' ${draft}_DpnII.txt > ${draft}.chrom.sizes
$juicer/scripts/juicer.sh -g ${draft} -s DpnII -z ${draft}.fa -y ./${draft}_DpnII.txt -S early -p ./${draft}.chrom.sizes -t 16
mkdir 3d_r3
cd 3d_r3
$d/run-asm-pipeline.sh -r 3 ./${draft}.fa ../aligned/merged_nodups.txt
echo "step1 done"
$d/run-asm-pipeline-post-review.sh -r ${draft}.rawchrom.review.assembly ../${draft}.fa ../aligned/merged_nodups.txt 
echo "3ddna done"