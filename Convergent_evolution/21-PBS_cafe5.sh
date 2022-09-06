#PBS -N cafe5
#PBS -l nodes=1:ppn=20
#PBS -q comput

PBS_O_WORKDIR=/public/home/WHC_sunyx/zx/Sobv_genome/CAFE5-master

cd $PBS_O_WORKDIR

./bin/cafe5  -i gene_families_filter2.txt -t tree.txt -p -k 3  -c 40 -o k3_re 
echo "k3 finished"
 
#./bin/cafe5  -i gene_families_filter2.txt -t tree.txt -p -k 4  -c 8 -o k4_re 
echo "k4 finished"

./bin/cafe5  -i gene_families_filter2.txt -t tree.txt -p -k 5  -c 40 -o k5_re 
echo "k4 finished"


./bin/cafe5  -i gene_families_filter2.txt -t tree.txt -p -k 6 -c 40 -o k6_re 
echo "k4 finished"


#############

#用orthofinder2的结果文件Orthogroups.GeneCount.tsv转换成gene_families.txt文件
#awk -v OFS="\t" '{$NF=null;print $1,$0}' Orthogroups.GeneCount.tsv |sed -E -e 's/Orthogroup/desc/' -e 's/_[^\t]+//g' >gene_families.txt
#
#如果是MCL获得的基因家族数据。cafe5的tutorials目录下有脚本mcl2rawcafe.py用于把mcl的输出转化成cafe5的输入。
#python mcl2rawcafe.py -i dump.blast_output.mci.I30 -o gene_families.txt -sp "ENSG00 ENSPTR ENSPPY ENSPAN ENSNLE ENSMMU ENSCJA ENSRNO ENSMUS ENSFCA ENSECA ENSBTA"。
#
#剔除不同物种间的基因拷贝数变异特别大的基因家族。
#否则可能无法估计参数，运行中断报错Failed to initialize any reasonable values
#
#cafe5的tutorial目录下脚本clade_and_size_filter.py可以筛除一个或以上物种有超过100个基因拷贝的基因家族。
#python clade_and_size_filter.py -s -i gene_families.txt -o gene_familie_filter.txt
#脚本运行失败。
#
#也可以用awk命令来筛除3-13列中基因家族数量大于等于100的行，可以根据自己的数据更改。
#awk 'NR==1 || $3<100 && $4<100 && $5<100 && $6<100 && $7<100 && $8<100 && $9<100 && $10<100 && $11<100 && $12<100 && $13<100  && $14<100 && $15<100 && $16<100 && $17<100  && $18<100 && $19<100 && $20<100  && $21<100  && $22<100 {print $0}' gene_families.txt >gene_families_filter.txt

#awk 'NR==1 || $3>0 && $4>0 && $5>0 && $6>0 && $7>0 && $8>0 && $9>0 && $10>0 && $11>0 && $12>0 && $13>0  && $14>0 && $15>0 && $16>0 && $17>0  && $18>0 && $19>0 && $20>0  && $21>0  && $22>0 {print $0}' gene_families_filter.txt >gene_families_filter2.txt

#用paml的Figtree转换成tree.txt文件（newick格式）
#grep -E -v "NEXUS|BEGIN|END" FigTree.tre|sed -E -e "s/\[[^]]*\]//g" -e "s/[ \t]//g" -e "/^$/d" -e "s/UTREE/tree tree/" >tree.txt

#对特定物种扩张和收缩基因的提取
#
#cat Gamma_change.tab |cut -f1,2|grep  -v "-"|sed '/0$/d'>Sobv.expanded #提取Gamma_change.tab第15列代表物种Sobv的扩张的orthogroupsID
#cat Gamma_change.tab |cut -f1,2|grep "-" >Sobv.contracted  #提取Gamma_change.tab第15列代表物种Sobv的收缩的orthogroupsID
#cat Gamma_family_results.txt |grep "y"|cut -f1 >p0.05.significant #提取显著扩张或收缩的orthogroupsID
#grep -f p0.05.significant Sobv.expanded |cut -f1>Sobv.expanded.significant #提取显著扩张的Sobv物种的orthogroupsID
#grep -f p0.05.significant Sobv.contracted |cut -f1 >Sobv.contracted.significant #提取显著收缩的Sobv物种的orthogroupsID
#
#grep -f Sobv.contracted.significant Orthogroups.txt |sed -E -e "s/: [^b]+bv/ bv/g" -e "s/ [^b]+//g" >Sobv.contracted.significant.ortho #提取显著收缩的基因
#grep -f Sobv.expanded.significant Orthogroups.txt |sed -E -e "s/: [^b]+bv/ bv/g" -e "s/ [^b]+//g" >Sobv.expanded.significant.ortho #提取显著扩张的基因
#
#cat Sobv.expanded.significant.ortho |sed "s/ /\n/g"|grep "bv" |sort -k 1.3n |uniq >Sobv.expanded.significant.genes
#cat Sobv.contracted.significant.ortho |sed "s/ /\n/g"|grep "bv" |sort -k 1.3n |uniq >Sobv.contracted.significant.genes
#
#seqkit grep -f Sobv.expanded.significant.genes Sobv.pep.fa >Sobv.expanded.significant.pep.fa
#seqkit grep -f Sobv.contracted.significant.genes Sobv.pep.fa >Sobv.contracted.significant.pep.fa
#提取出指定物种的显著扩张和收缩的蛋白序列之后，就可以拿去做GO注释和基因富集分析。
#
#3.4. 把每个节点收缩扩张的基因数量画在树上
#
#python python_scripts/cafetutorial_draw_tree.py -i reports/summary_run1_node.txt
# -t '((((cat:68.7105,horse:68.7105):4.56678,cow:73.2773):20.7227,(((((chimp:4.444
#17,human:4.44417):6.68268,orang:11.1268):2.28586,gibbon:13.4127):7.21153,(macaque
#:4.56724,baboon:4.56724):16.057):16.0607,marmoset:36.6849):57.3151):38.738,(rat:3
#6.3024,mouse:36.3024):96.4356)' -d '((((cat<0>,horse<2>)<1>,cow<4>)<3>,(((((chimp
#<6>,human<8>)<7>,orang<10>)<9>,gibbon<12>)<11>,(macaque<14>,baboon<16>)<15>)<13>,
#marmoset<18>)<17>)<5>,(rat<20>,mouse<22>)<21>)<19>' -o reports/summary_run1_tree_
#rapid.png -y Rapid
#
##########


