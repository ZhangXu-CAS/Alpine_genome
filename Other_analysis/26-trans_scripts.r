
####差异分析
df <- read.table("gene_count_matrix.csv", header = T,
                 sep = ",",
                 stringsAsFactors = F)
expr_df <- df[,c(1,13:21)]

coldata <- data.frame(condition=factor(c( "contral","contral","contral",
                                          "contral","contral","contral",
                                          "treat","treat","treat",
                                          "treat","treat","treat" )))

##DESeq2差异分析
#BiocManager::install("DESeq2")

library(DESeq2)

dds <- DESeqDataSetFromMatrix(df,
                       colData = coldata,
                       design = ~ condition,
                       tidy = TRUE)
dds <- DESeq(dds)
res <- results(dds)
dds_fliter<- na.omit(dds)

##PCA比较下组内差异是不是比组间差异小

rld <- rlog(dds)
plotPCA(rld)
##PCA
#dds <- estimateSizeFactors(dds)
#raw <- SummarizedExperiment(counts(dds, normalized=FALSE),
#                            colData=colData(dds))
#nor <- SummarizedExperiment(counts(dds, normalized=TRUE),
#                            colData=colData(dds))
#vsd <- vst(dds)
#rld <- rlog(dds)
#pdf("PCA.pdf")
#plotPCA( DESeqTransform(raw), intgroup=c("condition", "type") )
#plotPCA( DESeqTransform(nor), intgroup=c("condition", "type") )
#plotPCA(vsd, intgroup=c("condition", "type"))
#plotPCA(rld, intgroup=c("condition", "type"))
#dev.off()

library(ggplot2)


###绘制火山图#######
#####
DEG_data <- data.frame(res)
DEG_data$logP <- -log10(DEG_data$padj) # 对差异基因矫正后p-value进行log10()转换
dim(DEG_data)
DEG_data_filter<- na.omit(DEG_data[,-10:-11])
## [1] 5268    7
#将基因分为三类：not-siginficant，up，dowm
#将adj.P.value小于0.05，logFC大于1的基因设置为显著上调基因
#将adj.P.value小于0.05，logFC小于-1的基因设置为显著上调基因
DEG_data$Group <- "not-siginficant"
DEG_data$Group[which((DEG_data$padj < 0.05) & DEG_data$log2FoldChange > 1)] = "up-regulated"
DEG_data$Group[which((DEG_data$padj < 0.05) & DEG_data$log2FoldChange < -1)] = "down-regulated"
table(DEG_data$Group)
## 

DEG_data <- DEG_data[order(DEG_data$padj),]#对差异表达基因调整后的p值进行排序
#火山图中添加点(数据构建)
up_label <- head(DEG_data[DEG_data$Group == "up-regulated",],1)
down_label <- head(DEG_data[DEG_data$Group == "down-regulated",],1)
deg_label_gene <- data.frame(gene = c(rownames(up_label),rownames(down_label)),
                             label = c(rownames(up_label),rownames(down_label)))
DEG_data$gene <- rownames(DEG_data)
DEG_data <- merge(DEG_data,deg_label_gene,by = 'gene',all = T)
#不添加label
library(ggpubr)
ggscatter(DEG_data,x = "log2FoldChange",y = "logP",
          color = "Group",
          palette = c("blue","gray","red"),
          repel = T,
          ylab = "-log10(Padj)",
          size = 1) + 
  scale_y_continuous(limits = c(0,40))+
  scale_x_continuous(limits = c(-10,10))+
  geom_hline(yintercept = 1.3,linetype = "dashed",size=0.5)+
  geom_vline(xintercept = c(-1,1),linetype = "dashed",size=0.5)


 ####筛选差异表达基因
 
 
res_sub <- data.frame(subset(res, abs(log2FoldChange) >1& padj <0.05))
gene_up <- row.names(res_sub[res_sub$log2FoldChange > 0, ])
length(gene_up) # 1024
gene_down <- row.names(res_sub[res_sub$log2FoldChange < 0, ])
length(gene_down) # 702

deg_gene<- c(gene_up,gene_down)
write.csv(gene_up,file = "Rh_DESeq_up.csv",row.names = F)
write.csv(gene_down,file = "Rh_DESeq_down.csv",row.names = F)
res_diff_data <- merge(as.data.frame(res_sub),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
write.csv(res_diff_data,file = "831_DESeq_diff_genes.csv",row.names = F)

res_data <- merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
write.csv(res_data,file = "8.31_DESeq_all_genes.csv",row.names = F)


###聚类热图
library(genefilter)
library(pheatmap)
library(GiNA)
rld <- rlog(dds)
#rld <- rlogTransformation(dds_out,blind = F)
#write.csv(assay(rld),file="mm.DESeq2.pseudo.counts.csv")

#topVarGene <- head(order(rowVars(assay(rld)),decreasing = TRUE),2390)
mat  <- assay(rld)[rownames(res_sub),]
mat <- mat - rowMeans(mat) #减去一个平均值，让数值更加集中。第二个图
dim(mat)
anno <- as.data.frame(colData(rld)[,c("condition","sizeFactor")])
pheatmap(mat,color = colorRampPalette(c("navy", "white", "firebrick3"))(100), 
         #color = colorRampPalette(rev(brewer.pal(n = 7, name =  "RdYlBu")))(100),
         show_rownames=0,
         clustering_distance_rows="euclidean",
         scale = "row"
         )

####使用clusterProfiler做GO/KEGG富集分析

install.packages("https://github.com/xuzhougeng/org.Osativa.eg.db/releases/download/v0.01/org.Osativa.eg.db.tar.gz", 
                 repos = NULL, 
                 type="source")

library(clusterProfiler)
library(org.Osativa.eg.db)
#> query(hub, "oryza sativa")
# AnnotationHub with 6 records
# snapshotDate(): 2021-05-18
# $dataprovider: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/, WikiPathways, NCBI,DBCLS, Inparanoid8
# $species: Oryza sativa, Oryza sativa_subsp._japonica, Oryza sativa_Japonica_Group, Oryza sa...
# $rdataclass: OrgDb, Tibble, SQLiteFile, Inparanoid8Db
# additional mcols(): taxonomyid, genome, description, coordinate_1_based, maintainer,
#   rdatadateadded, preparerclass, tags, rdatapath, sourceurl, sourcetype 
# retrieve records with, e.g., 'object[["AH10561"]]' 

#title                                               
#AH10561 | hom.Oryza_sativa.inp8.sqlite                        
#AH91623 | MeSHDb for Oryza sativa Japonica Group (Rice, v001) 
#AH91809 | wikipathways_Oryza_sativa_metabolites.rda           
#AH94059 | org.Oryza_sativa_(japonica_cultivar-group).eg.sqlite
#AH94060 | org.Oryza_sativa_Japonica_Group.eg.sqlite           
#AH94061 | org.Oryza_sativa_subsp._japonica.eg.sqlite
#org <- hub[['AH94060']]
#gene_up2<- bitr(gene_up, fromType, toType, OrgDb)

org <- org.Osativa.eg.db
ego_up <- enrichGO(gene_up, pvalueCutoff = 0.05, 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.2, 
                   OrgDb = org,
                   keyType = "GID",
                   ont="all")
 dotplot(ego_up)

ego_down <- enrichGO(gene_down, pvalueCutoff = 0.05, 
                     OrgDb = org,
                     keyType = "GID",
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.2, 
                     ont="all")
dotplot(ego_down)
ego_both <- enrichGO(deg_gene,
                   OrgDb = org,
                   keyType = "GID",
                   pAdjustMethod = "none",
                   ont="all")
dotplot(ego_both)





