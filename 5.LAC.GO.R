library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)

pvalueFilter=0.05      
qvalueFilter=0.05       

colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}
ontology.col=c("#00AFBB", "#E7B800", "#90EE90")

rt=read.table("diffGene.txt", header=T, sep="\t", check.names=F)  

genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"] 

kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
write.table(GO, file="GO.txt", sep="\t", quote=F, row.names = F)

showNum=6
if(nrow(GO)<18){
	showNum=nrow(GO)
}

pdf(file="barplot.pdf", width=9, height=12)
bar=barplot(kk, drop=TRUE, showCategory=showNum, label_format=30, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()
		