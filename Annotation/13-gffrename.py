#!/usr/bin/env python3
#ref. https://blog.csdn.net/u012110870/article/details/82500684, 如何对基因组序列进行注释
##from 徐洲更hoptop
import re
import sys

if len(sys.argv) < 3:
    sys.exit()

gff = open(sys.argv[1])
prf = sys.argv[2]

count = 0
mRNA  = 0
cds   = 0
exon  = 0

print("##gff-version 3.2.1")
for line in gff:
    if not line.startswith("\n"):
        records = line.split("\t")
        records[1] = "."
    if re.search(r"\tgene\t", line):
        count = count + 1
        mRNA  = 0
        gene_id = prf + str(count).zfill(5)
        records[8] = "ID={}".format(gene_id)
    elif re.search(r"\tmRNA\t", line):
        cds   = 0
        exon  = 0
        mRNA  = mRNA + 1
        mRNA_id    = gene_id + "." + str(mRNA)
        records[8] = "ID={};Parent={}".format(mRNA_id, gene_id)
    elif re.search(r"\texon\t", line):
        exon     = exon + 1
        exon_id  = mRNA_id + "_exon_" + str(exon)
        records[8] = "ID={};Parent={}".format(exon_id, mRNA_id)
    elif re.search(r"\tCDS\t", line):
        cds     = cds + 1
        cds_id  = mRNA_id + "_cds_" + str(cds)
        records[8] = "ID={};Parent={}".format(cds_id, mRNA_id)
    else:
        continue

    print("\t".join(records))

gff.close()
