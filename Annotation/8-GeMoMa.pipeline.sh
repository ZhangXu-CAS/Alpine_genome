#! ~/bin/bash

java -jar GeMoMa/GeMoMa-1.7.1.jar CLI GeMoMaPipeline threads=24 \
outdir=chr.masked+trans GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=chr.masked.fa\
	 s=own i=Athaliana g=Athaliana_447_TAIR10.fa\
	 a=Athaliana_167_TAIR10.gene_exons.gff3\
	 s=own i=Osativa g=Osativa_323_v7.0.fa\
	 a=Osativa_323_v7.0.gene_exons.gff3\
	 s=own i=Stuberosum g=Stuberosum_448_v4.03.fa\
	 a=Stuberosum_448_v4.03.gene_exons.gff3\
	 s=own i=Helianthus g=GCF_002127325.2_HanXRQr2.0-SUNRISE_genomic.fna\
	 a=GCF_002127325.2_HanXRQr2.0-SUNRISE_genomic.gff\
	 s=own i=Lactuca  g=GCF_002870075.1_Lsat_Salinas_v7_genomic.fna\
	 a=GCF_002870075.1_Lsat_Salinas_v7_genomic.gff\
	 s=own i=Cynara  g=GCF_001531365.1_CcrdV1_genomic.fna\
	 a=GCF_001531365.1_CcrdV1_genomic.gff\
	 s=own i=Mikania  g=GCA_009363875.1_ASM936387v1_genomic.fna\
	 a=GCA_009363875.1_ASM936387v1_genomic.gff\
	 r=MAPPED ERE.m=/home/wangh/zx/Saussurea_genome/annotation/trans_evi/merged.bam

echo 'geoma done'

