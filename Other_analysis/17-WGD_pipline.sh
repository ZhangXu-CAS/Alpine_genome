#!/bin/bash

####Panparalogs
for i in 'Cynara' 'Helianthus' 'Sobv'

do 
{
wgd mcl --cds --mcl -s ${i}_cds.fa -o ${i}.paralogs -n 12
wgd ksd -o $i.wgd  $i.paralogs/*.mcl ${i}_cds.fa -n 12
wgd mix $i.wgd/*ks.tsv -o $i.mix -n 1 5
} &
 done 
 wait
 
#wgd mcl --cds --mcl -s Sobv_v1_cds.fa -o Sobv_v1.paralogs -n 16
#wgd ksd -o Sobv_v1.wgd  Sobv_v1.paralogs/*.mcl Sobv_v1_cds.fa  -n 16
#wgd mix Sobv.wgd/*ks.tsv -o Sobv.mix -n 1 5

 wgd syn --feature mRNA  --gene_attribute ID -o Sobv.syn -ks Sobv_cds.fa.ks.tsv Sobv.gff3  Sobv_cds.fa.blast.tsv.mcl


###Orthologs

for i in 'Cynara' 'Helianthus' 

do 
{
#mkdir Sobv_$i.out 
wgd mcl --cds --one_v_one -s ${i}_cds.fa,Sobv_v1_cds.fa -id $i,Sobv -o Sobv_v1.$i.out  -n 10
wgd ksd -o $i.Sobv_v1.orthologs Sobv_v1.$i.out/*.ovo.tsv ${i}_cds.fa Sobv_v1_cds.fa  -n 10

 } &
 done 
 wait
 
#wgd viz -ks ./Ks_tsvs -i 
#wgd viz -ks ./Ks_tsvs -a 0.5,0.5,0.5,0.5,0.5 -c Salmon,Gold,Green,RoyalBlue,Yellow -b 100 -ht barstacked
#wgd viz -ks Ks_tsvs  -i  -l Sunflower,Sunflower-Saussurea,Artichoke,Artichoke-Saussrea,Saussurea



