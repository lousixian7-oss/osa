
library(limma)
library(ggpubr)

LACCluFile="LACcluster.txt"       
geneCluFile="geneCluster.txt"      
scoreFile="LACscore.txt"           

LACClu=read.table(LACCluFile, header=T, sep="\t", check.names=F, row.names=1)
geneClu=read.table(geneCluFile, header=T, sep="\t", check.names=F, row.names=1)
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)

twoCluster=cbind(LACClu, geneClu)
sameSample=intersect(row.names(twoCluster), row.names(score))
data=cbind(score[sameSample,,drop=F], twoCluster[sameSample,c("LACcluster","geneCluster"),drop=F])

data$LACcluster=factor(data$LACcluster, levels=levels(factor(data$LACcluster)))
group=levels(factor(data$LACcluster))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data$LACcluster)))]
	
boxplot=ggboxplot(data, x="LACcluster", y="LACscore", color="LACcluster",
			      xlab="LACcluster",
			      ylab="LACscore",
			      legend.title="LACcluster",
			      palette=bioCol,
			      add = "jitter")+ 
	stat_compare_means(comparisons = my_comparisons)
	
pdf(file="LACcluster.pdf", width=5, height=4.5)
print(boxplot)
dev.off()

data$geneCluster=factor(data$geneCluster, levels=levels(factor(data$geneCluster)))
group=levels(factor(data$geneCluster))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data$geneCluster)))]
	
boxplot=ggboxplot(data, x="geneCluster", y="LACscore", color="geneCluster",
			      xlab="geneCluster",
			      ylab="LACscore",
			      legend.title="geneCluster",
			      palette=bioCol,
			      add = "jitter")+ 
	stat_compare_means(comparisons = my_comparisons)
	
pdf(file="geneCluster.pdf", width=5, height=4.5)
print(boxplot)
dev.off()
