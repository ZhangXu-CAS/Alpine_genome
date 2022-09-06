#!/usr/bin/bash

## Pairwise Sequentially Markovian Coalescent (PSMC) model
## Ref. web, https://github.com/lh3/psmc
## Write by Xu Zhang, 2020/2/29 10:00am



REF=/home/wangh/zx/Saussurea_genome/depth/chr.fa

BAM=/home/wangh/zx/Saussurea_genome/depth/chr.sort.mark.bam
 samtools mpileup -C50 -uf $REF $BAM | bcftools view -c - \
      | vcfutils.pl vcf2fq -d 20 -D 120 | gzip > chr_psmc.fq.gz
##Here option -d sets and minimum read depth and -D sets the maximum. It is recommended to set -d to a third of the average depth and -D to twice. 
./utils/fq2psmcfa -q20 chr_psmc.fq.gz > chr_psmc.psmcfa 
##transforms the consensus sequence into a fasta-like format where the i-th character in the output sequence.
 
./psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o chr_psmc.psmc chr_psmc.psmcfa

#The `psmc' program infers the scaled mutation rate, the recombination rate and the free population size parameters.


#./utils/psmc2history.pl chr_psmc.psmc | ./utils/history2ms.pl > ms-cmd.sh

#`psmc2history.pl' combined with `history2ms.pl' to generate the ms command line that simulates the history inferred by PSMC, or visualize the result with `psmc_plot.pl'.

./utils/psmc_plot.pl chr_psmc  -u xx -g 5 chr_psmc.psmc

done
echo "done"


#To perform bootstrapping:

./utils/splitfa chr_psmc.psmcfa > split.psmcfa
seq 100 | xargs -i echo ./psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o round-{}.psmc split.psmcfa | sh
cat chr_psmc.psmc round-*.psmc > chr_psmc.combined.psmc
./utils/psmc_plot.pl -u 6.59e-09 -g 6 combinedchr_psmc.combined.psmc 

#Correcting for low coverage

#	python $WD psmc_plot.pl -M "sample1=0.1,sample2=0.2" prefix sample1.psmc sample2.psmc 

##This says that sample1 has 10% false negative rate (FNR) on hets and sample2has 20%.