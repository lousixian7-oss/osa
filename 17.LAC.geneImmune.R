
library(limma)
library(pheatmap)
library(reshape2)
library(ggpubr)

gene="TIMM50"                       #PNPO,RMND1
clusterFile="LACcluster.txt"        
ssgseaFile="ssGSEA.result.txt"      

cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

ssgsea=read.table(ssgseaFile, header=T, sep="\t", check.names=F, row.names=1)
ssgsea=t(ssgsea)

sameSample=intersect(row.names(cluster), row.names(ssgsea))
cluster=cluster[sameSample,gene,drop=F]
ssgsea=ssgsea[sameSample,,drop=F]
rt=cbind(ssgsea, cluster)
rt[,gene]=ifelse(rt[,gene]>median(rt[,gene]), "High", "Low")
rt[,gene]=factor(rt[,gene], levels=c("Low", "High"))

data=melt(rt, id.vars=c(gene))
colnames(data)=c("Gene", "Immune", "Fraction")

bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data[,"Gene"])))]
p=ggboxplot(data, x="Immune", y="Fraction", color="Gene",
     xlab="",
     ylab="Immune infiltration",
     legend.title=gene,
     palette=bioCol)
p=p+rotate_x_text(50)
#????ͼƬ?ļ?
pdf(file="boxplotTIMM50.pdf", width=8, height=6.5)
p+stat_compare_means(aes(group=Gene),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")
dev.off()

