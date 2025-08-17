
library(ggalluvial)
library(ggplot2)
library(dplyr)

LACCluFile="LACcluster.txt"      
geneCluFile="geneCluster.txt"     
scoreFile="LACscore.txt"          

LACClu=read.table(LACCluFile, header=T, sep="\t", check.names=F, row.names=1)
geneClu=read.table(geneCluFile, header=T, sep="\t", check.names=F, row.names=1)
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)

twoCluster=cbind(LACClu, geneClu)
sameSample=intersect(row.names(twoCluster), row.names(score))
rt=cbind(twoCluster[sameSample,c("LACcluster","geneCluster"),drop=F], score[sameSample,,drop=F])
rt$LACscore=ifelse(rt$LACscore>0, "High", "Low")

corLodes=to_lodes_form(rt, axes = 1:ncol(rt), id = "Cohort")

pdf(file="ggalluvial.pdf", width=6, height=6)
mycol=rep(c("#0066FF","#FF0000","#FF9900","#029149","#6E568C","#E0367A","#D8D155","#223D6C","#D20A13","#431A3D","#91612D","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#64495D","#7CC767"),15)
ggplot(corLodes, aes(x = x, stratum = stratum, alluvium = Cohort,fill = stratum, label = stratum)) +
  	 scale_x_discrete(expand = c(0, 0)) +  
  	 geom_flow(width = 2/10, aes.flow = "forward") + 
	 geom_stratum(alpha = .9,width = 2/10) +
	 scale_fill_manual(values = mycol) +
	 geom_text(stat = "stratum", size = 3,color="black") +
	 xlab("") + ylab("") + theme_bw() + 
	 theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank()) + #ȥ????????
	 theme(panel.grid =element_blank()) + 
	 theme(panel.border = element_blank()) + 
	 ggtitle("") + guides(fill = FALSE)                            
dev.off()
