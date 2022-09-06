#从gff3文件获取基因的位置
grep '[[:blank:]]mRNA[[:blank:]]' S.obv.gff3 | cut -f 1,4,5 | awk '{print $1"\t"$2"\t"$3}'|grep chr > genes.bed
#从基因组位置文件获取滑动窗
bedtools makewindows -g chr.fa.fai -w 1000000 >chr.windows 
#获取基因覆盖率
bedtools coverage -a chr.windows -b genes.bed | cut -f 1-4 > genes_num.txt

#获取GC含量
bedtools nuc -fi chr.fa -bed chr.windows |cut -f 1-3,5|sed '1d' >GC_content_num.txt
###Copia类型的转座子
cat Spohua.fa.pass.list|grep Copia |cut -f1|tr : "\t"|sed 's/\.\./\t/g' >chr.Copia.bed
###Gypsy类型的转座子
cat Spohua.fa.pass.list|grep Gypsy |cut -f1|tr : "\t"|sed 's/\.\./\t/g' >chr.Gypsy.bed
##获取覆盖率

bedtools coverage -a chr.windows -b chr.Copia.bed |cut -f 1-4 >Copia_num.txt
bedtools coverage -a chr.windows -b chr.Gypsy.bed |cut -f 1-4 >Gypsy_num.txt

###完整LTR的百分比（每个windowsize的完整LTR的百分比）
cat Spohua.fa.out.LAI|awk '{print $1"\t"$2"\t"$3"\t"$4*100}'|sed '1,2d' >full_LTR_num.txt

###重复序列
grep '[[:blank:]]*_retrotransposon[[:blank:]]' chr.fa.mod.EDTA.intact.gff3 | cut -f 1,4,5 | awk '{print $1"\t"$2"\t"$3}'|grep chr > LTR.bed
bedtools coverage -a chr.windows -b LTR.bed | cut -f 1-4 > LTR_num.txt

grep '[[:blank:]]Gypsy_LTR_retrotransposon[[:blank:]]' chr.fa.mod.EDTA.intact.gff3 | cut -f 1,4,5 | awk '{print $1"\t"$2"\t"$3}'|grep chr > Gypsy.bed
bedtools coverage -a chr.windows -b Gypsy.bed | cut -f 1-4 > Gypsy_num.txt

grep '[[:blank:]]Copia_LTR_retrotransposon[[:blank:]]' chr.fa.mod.EDTA.intact.gff3 | cut -f 1,4,5 | awk '{print $1"\t"$2"\t"$3}'|grep chr > Copia.bed
bedtools coverage -a chr.windows -b Copia.bed | cut -f 1-4 > Copia_num.txt

grep '[[:blank:]]LTR_retrotransposon[[:blank:]]' chr.fa.mod.EDTA.intact.gff3 | cut -f 1,4,5 | awk '{print $1"\t"$2"\t"$3}'|grep chr > LTR_retrotransposon.bed
bedtools coverage -a chr.windows -b LTR_retrotransposon.bed | cut -f 1-4 > Copia_num.txt