
library(limma)
library(ggplot2)

clusterFile="LACcluster.txt"     

rt=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
data=rt[,1:(ncol(rt)-1),drop=F]
LACcluster=as.vector(rt[,ncol(rt)])

data.pca=prcomp(data)
pcaPredict=predict(data.pca)
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], LACcluster=LACcluster)
PCA.mean=aggregate(PCA[,1:2], list(LACcluster=PCA$LACcluster), mean)

bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
LACCluCol=bioCol[1:length(levels(factor(LACcluster)))]

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
}
df_ell <- data.frame()
for(g in levels(factor(PCA$LACcluster))){
df_ell <- rbind(df_ell, cbind(as.data.frame(with(PCA[PCA$LACcluster==g,],
                  veganCovEllipse(cov.wt(cbind(PC1,PC2),
                  wt=rep(1/length(PC1),length(PC1)))$cov,
                  center=c(mean(PC1),mean(PC2))))), LACcluster=g))
}

pdf(file="PCA.pdf", height=5, width=6.5)
ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = LACcluster)) +
	scale_colour_manual(name="LACcluster", values =LACCluCol)+
    theme_bw()+
    theme(plot.margin=unit(rep(1.5,4),'lines'))+
    geom_path(data=df_ell, aes(x=PC1, y=PC2, colour=LACcluster), size=1, linetype=2)+
    annotate("text",x=PCA.mean$PC1, y=PCA.mean$PC2, label=PCA.mean$LACcluster, cex=7)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
