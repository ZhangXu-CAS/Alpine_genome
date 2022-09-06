#!~/bin/bash:
WD=/home/wangh/zx/Saussurea_genome/annotation/trans_evi
cd $WD

mkdir index
hisat2-build chr.fa index/chr
hisat2 --dta -p 12 -x index/chr -1 BX-L_clean_1.fq.gz -2 BX-L_clean_2.fq.gz| samtools sort -@ 10 > BX-L.bam & 
hisat2 --dta -p 12 -x index/chr -1 BX-F_clean_1.fq.gz -2 BX-F_clean_2.fq.gz| samtools sort -@ 10 > BX-F.bam &
hisat2 --dta -p 12 -x index/chr -1 BX-B_clean_1.fq.gz -2 BX-B_clean_2.fq.gz| samtools sort -@ 10 > BX-B.bam &
samtools merge -@ 20 merged.bam BX-L.bam BX-F.bam  BX-B.bam 
stringtie -p 20 -o HISAT2+StringTie.gtf merged.bam

####TransDecoder进行编码区预测
./util/gtf_genome_to_cdna_fasta.pl HISAT2+StringTie.gtf ../chr.fa > transcripts.fasta
./util/gtf_to_alignment_gff3.pl HISAT2+StringTie.gtf > transcripts.gff3
./TransDecoder.LongOrfs -t transcripts.fasta
./TransDecoder.Predict -t transcripts.fasta
./util/cdna_alignment_orf_to_genome_orf.pl \
     transcripts.fasta.transdecoder.gff3 \
     transcripts.gff3 \
     transcripts.fasta > transcripts.fasta.transdecoder.genome.gff3
echo 'TransDecoder done'

###最后结果transcripts.fasta.transdecoder.gff3用于提供给EvidenceModeler