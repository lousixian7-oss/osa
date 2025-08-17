
library(limma)
library(pheatmap)
library(reshape2)
library(ggpubr)

clusterFile="LACcluster.txt"         
ssgseaFile="ssGSEA.result.txt"       

cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
cluster=cluster[,-ncol(cluster)]

ssgsea=read.table(ssgseaFile, header=T, sep="\t", check.names=F, row.names=1)
ssgsea=t(ssgsea)

sameSample=intersect(row.names(cluster), row.names(ssgsea))
cluster=cluster[sameSample,,drop=F]
ssgsea=ssgsea[sameSample,,drop=F]

cor=cor(ssgsea, cluster, method="spearman")

#????????????ͼ
pdf(file="heatmap.pdf", width=7, height=4.5)
pheatmap(cor,
         color = colorRampPalette(c(rep("blue",1), "white", rep("red",1)))(100),
         cluster_cols =F,
         cluster_rows =T,
         display_numbers = T,
         show_colnames=T,
         show_rownames=T,
         angle_col =45,
         fontsize=6,
         fontsize_row=7,
         fontsize_col=6)
dev.off()
