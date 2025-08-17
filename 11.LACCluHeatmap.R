library(limma)
library(pheatmap)
library(reshape2)
library(ggpubr)

clusterFile="LACcluster.txt"      

rt=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
rt=rt[order(rt$LACcluster),]

data=t(rt[,1:(ncol(rt)-1),drop=F])
Type=rt[,ncol(rt),drop=F]

bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
ann_colors=list()
LACCluCol=bioCol[1:length(levels(factor(Type$LACcluster)))]
names(LACCluCol)=levels(factor(Type$LACcluster))
ann_colors[["LACcluster"]]=LACCluCol

pdf("heatmap.pdf", width=7, height=4.5)
pheatmap(data,
         annotation=Type,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("blue",2), "white", rep("red",2)))(100),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=6,
         fontsize_row=7,
         fontsize_col=6)
dev.off()

data=melt(rt, id.vars=c("LACcluster"))
colnames(data)=c("LACcluster", "Gene", "Expression")

p=ggboxplot(data, x="Gene", y="Expression", color = "LACcluster", 
	     ylab="Gene expression",
	     xlab="",
	     legend.title="LACcluster",
	     palette = LACCluCol,
	     width=1)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=LACcluster),
	      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
	      label = "p.signif")

pdf(file="boxplot.pdf", width=6, height=5)
print(p1)
dev.off()

