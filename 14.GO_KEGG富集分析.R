library(clusterProfiler)#GO富集分析、KEGG通路富集分析
library(org.Hs.eg.db)#基因注释数据库
library(enrichplot)
library(ggplot2)
library(GOplot)
#设置工作路径
setwd("工作路径")
#读入数据
genes_df <- read.csv("genes.csv")  

# 将基因符号转换为ENTREZ ID
entrezIDs <- bitr(genes_df$gene, fromType = "SYMBOL", 
               toType = c("ENTREZID", "SYMBOL"),
               OrgDb = org.Hs.eg.db) # `OrgDb`参数指定使用的人类基因组注释数据库

# 去掉前后空格
genes_df$gene <- trimws(genes_df$gene)

# 查看没匹配上的基因
library(dplyr)
mapped <- bitr(genes_df$gene, 
               fromType = "SYMBOL", 
               toType = c("ENTREZID", "SYMBOL"),
               OrgDb = org.Hs.eg.db)

unmapped <- setdiff(genes_df$gene, mapped$SYMBOL)
unmapped

# 可用别名查询
mapped_alias <- bitr(unmapped, 
                     fromType = "ALIAS", 
                     toType = c("ENTREZID", "SYMBOL"), 
                     OrgDb = org.Hs.eg.db)
mapped_symbol <- bitr(genes_df$gene,
                      fromType = "SYMBOL",
                      toType = c("ENTREZID", "SYMBOL"),
                      OrgDb = org.Hs.eg.db)

entrezIDs <- bind_rows(mapped_symbol, mapped_alias) %>%
  distinct()#使用entrezIDs 
gene<- entrezIDs$ENTREZID
##GO富集分析
go<- enrichGO(gene = gene,OrgDb = org.Hs.eg.db, pvalueCutoff =1, qvalueCutoff = 1,ont="all",readable =T)
write.table(go,file="GO.txt",sep="\t",quote=F,row.names = F) #

##可视化
##条形图
pdf(file="GO-柱状图.pdf",width = 8,height = 7)
##showCategory改变自己想要展示的条目数量
barplot(go, drop = TRUE, showCategory =3,label_format=100,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

##气泡图
pdf(file="GO-气泡图.pdf",width = 8,height = 7)
dotplot(go,showCategory = 6,label_format=100,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()


#kegg分析
kk <- enrichKEGG(gene = gene,keyType = "kegg",organism = "hsa", pvalueCutoff =1, qvalueCutoff =1, pAdjustMethod = "fdr")   
write.table(kk,file="KEGG.txt",sep="\t",quote=F,row.names = F)                         
kk@result <- kk@result[kk@result$Description != "", ]

##可视化
##条形图
pdf(file="KEGG-柱状图.pdf",width = 8,height = 6)
barplot(kk, drop = TRUE, showCategory = 10,label_format=100)
dev.off()

##气泡图
pdf(file="KEGG-气泡图.pdf",width = 8,height = 6)
dotplot(kk, showCategory = 8,label_format=100)
dev.off()
