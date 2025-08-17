
library(limma)
library(pheatmap)

expFile="interGeneExp.txt"       
clusterFile="geneCluster.txt"     

exp=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
exp=t(exp)

cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

sameSample=intersect(row.names(exp), row.names(cluster))
exp=exp[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
rt=cbind(exp, cluster)
rt=rt[order(rt$geneCluster),]

data=t(rt[,1:(ncol(rt)-1),drop=F])
Type=rt[,ncol(rt),drop=F]

bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
ann_colors=list()
geneCluCol=bioCol[1:length(levels(factor(Type$geneCluster)))]
names(geneCluCol)=levels(factor(Type$geneCluster))
ann_colors[["geneCluster"]]=geneCluCol

pdf(file="heatmap.pdf", width=8, height=6)
pheatmap(data,
         annotation_col =  Type,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("blue",2), "white", rep("red",2)))(100),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=6,
         fontsize_row=6,
         fontsize_col=6)
dev.off()

