

expFile="diffGeneExp.txt"     

data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)

pca=prcomp(data, scale=TRUE)
value=predict(pca)
LACscore=value[,1]+value[,2]
LACscore=as.data.frame(LACscore)
scoreOut=rbind(id=colnames(LACscore), LACscore)
write.table(scoreOut, file="LACscore.txt", sep="\t", quote=F, col.names=F)


