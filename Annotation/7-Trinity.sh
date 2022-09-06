
Trinity --genome_guided_bam merged.bam --max_memory 100G --genome_guided_max_intron 10000 --CPU 80 
Trinity --seqType fq --CPU 40 --max_memory 64G --left BX-B_clean_1.fq.gz,BX-L_clean_1.fq.gz,BX-F_clean_1.fq.gz --right BX-B_clean_2.fq.gz,BX-L_clean_2.fq.gz,BX-F_clean_2.fq.gz --no_bowtie 
echo 'trinity done'

