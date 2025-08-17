
library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)

expFile="ad.normalize.txt"         
gmtFile="immune.gmt"            
clusterFile="LACcluster.txt"    

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())

gsvaP = ssgseaParam(exprData = data, geneSets = geneSets, assay = NA_character_,annotation = NULL, minSize = 1, maxSize = Inf,alpha = 0.25, normalize = TRUE)
ssgseaScore= gsva(gsvaP)

normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)

ssgseaOut=rbind(id=colnames(ssgseaScore), ssgseaScore)
write.table(ssgseaOut,file="ssGSEA.result.txt",sep="\t",quote=F,col.names=F)

cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

ssgseaScore=t(ssgseaScore)
sameSample=intersect(row.names(ssgseaScore), row.names(cluster))
ssgseaScore=ssgseaScore[sameSample,,drop=F]
cluster=cluster[sameSample,"LACcluster",drop=F]
scoreCluster=cbind(ssgseaScore, cluster)

data=melt(scoreCluster, id.vars=c("LACcluster"))
colnames(data)=c("LACcluster", "Immune", "Fraction")

bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data[,"LACcluster"])))]
p=ggboxplot(data, x="Immune", y="Fraction", color="LACcluster",
            xlab="",
            ylab="Immune infiltration",
            legend.title="LACcluster",
            palette=bioCol)
p=p+rotate_x_text(50)

pdf(file="boxplot.pdf", width=8, height=6)
p+stat_compare_means(aes(group=LACcluster),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")
dev.off()
