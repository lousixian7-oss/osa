
library(rms)
library(rmda)

inputFile="ad.rfGeneExp.txt"      

data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
rt=cbind(as.data.frame(data), Type=group)
paste(colnames(data), collapse="+")

ddist=datadist(rt)
options(datadist="ddist")
ͼ
lrmModel=lrm(Type~ TIMM50+PNPO+RMND1, data=rt, x=T, y=T)
nomo=nomogram(lrmModel, fun=plogis,
	fun.at=c(0.0001,0.1,0.3,0.5,0.7,0.9,0.99),
	lp=F, funlabel="Risk of Disease")

pdf("Nom.pdf", width=8, height=6)
plot(nomo)
dev.off()

cali=calibrate(lrmModel, method="boot", B=1000)
pdf("Calibration.pdf", width=6, height=6)
plot(cali,
	xlab="Predicted probability",
	ylab="Actual probability", sub=F)
dev.off()

rt$Type=ifelse(rt$Type=="n", 0, 1)
dc=decision_curve(Type ~ TIMM50+PNPO+RMND1, data=rt, 
	family = binomial(link ='logit'),
	thresholds= seq(0,1,by = 0.01),
	confidence.intervals = 0.95)

pdf(file="DCA.pdf", width=6, height=6)
plot_decision_curve(dc,
	curve.names="LAC genes",
	xlab="Threshold probability",
	cost.benefit.axis=T,
	col="red",
	confidence.intervals=FALSE,
	standardize=FALSE)
dev.off()

pdf(file="clinical_impact.pdf", width=6, height=6)
plot_clinical_impact(dc,
	confidence.intervals=T,
	col = c("red", "blue"))
dev.off()

