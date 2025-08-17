
library(limma)
library(pheatmap)
library(reshape2)
library(ggpubr)

logFCfilter=0.3                
adj.P.Val.Filter=0.05      
inputFile="ad.LACGeneExp.txt"      
setwd("C:/Users/23634/Desktop/论文代码")

rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

Type <- sub("^.*_", "", colnames(data))

design <- model.matrix(~0+factor(Type))
colnames(design) <- c("n","o")
fit <- lmFit(data,design)
cont.matrix<-makeContrasts(o-n,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

allDiff=topTable(fit2,adjust='fdr',number=200000)
allDiffOut=rbind(id=colnames(allDiff),allDiff)
write.table(allDiffOut, file="all.txt", sep="\t", quote=F, col.names=F)

diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & adj.P.Val < adj.P.Val.Filter )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut, file="diff.txt", sep="\t", quote=F, col.names=F)

diffGeneExp=data[row.names(diffSig),]
diffGeneExpOut=rbind(id=paste0(colnames(diffGeneExp)), diffGeneExp)
write.table(diffGeneExpOut, file="diffGeneExp.txt", sep="\t", quote=F, col.names=F)

sigVec=c()
sigGeneVec=c()
for(i in row.names(diffSig)){
  adj.P.Val=diffSig[i,]$adj.P.Val
	Sig=ifelse(adj.P.Val<0.001,"***",ifelse(adj.P.Val<0.01,"**",ifelse(adj.P.Val<0.05,"*","")))
	if(adj.P.Val<0.05){
	sigVec=c(sigVec, paste0(i, Sig))
	sigGeneVec=c(sigGeneVec, i)}
}

row.names(diffGeneExp)=sigVec

names(Type)=colnames(diffGeneExp)
Type=as.data.frame(Type)
pdf("heatmap.pdf", width=7.5, height=4.7)
pheatmap(diffGeneExp,
         annotation=Type,
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

exp=as.data.frame(t(diffGeneExp))
exp=cbind(exp, Type=Type)
data=melt(exp, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")

p=ggboxplot(data, x="Gene", y="Expression", color = "Type", 
	     xlab="",
	     ylab="Gene expression",
	     legend.title="Type",
	     palette = c("blue", "red"),
	     width=1)
p=p+rotate_x_text(60)

pdf(file="boxplot.pdf", width=7.5, height=5)
print(p)
dev.off()


