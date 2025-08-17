
library(limma)
library(pheatmap)
library(reshape2)
library(ggpubr)


expFile="diffGeneExp.txt"        
clusterFile="geneCluster.txt"      

exp=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
exp=t(exp)

cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

sameSample=intersect(row.names(exp), row.names(cluster))
exp=exp[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
rt=cbind(exp, cluster)
rt=rt[order(rt$geneCluster),]

data=melt(rt, id.vars=c("geneCluster"))
colnames(data)=c("geneCluster", "Gene", "Expression")

bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data[,"geneCluster"])))]
p=ggboxplot(data, x="Gene", y="Expression", color = "geneCluster", 
	     xlab="",
	     ylab="Gene expression",
	     legend.title="geneCluster",
	     palette = bioCol,
	     width=1)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=geneCluster),
	      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
	      label = "p.signif")

pdf(file="boxplot.pdf", width=7, height=5)
print(p1)
dev.off()
