
library(pROC)    
library("ggplot2")

rt=read.table("ad.rfGeneExp.txt", header = TRUE, sep = "\t",
              check.names = FALSE, row.names = 1)
rt=t(rt)
suffix <- ifelse(grepl("_", rownames(rt)),
                 sub("^.*_", "", rownames(rt)),
                 NA_character_)   
rt <- data.frame(Group = suffix,  rt, check.names = FALSE) 

y=colnames(rt)[1]

bioCol=c("red","blue","green","yellow")
if(ncol(rt)>4){
  bioCol=rainbow(ncol(rt))}

pdf("ROC(多基因).pdf",width=5,height=5)

roc1=roc(rt[,y], as.vector(rt[,2]))
aucText=c( paste0(colnames(rt)[2],", AUC=",sprintf("%0.3f",auc(roc1))) )
plot(roc1, col=bioCol[1])

for(i in 3:ncol(rt)){
  roc1=roc(rt[,y], as.vector(rt[,i]))
  lines(roc1, col=bioCol[i-1])
  aucText=c(aucText, paste0(colnames(rt)[i],", AUC=",sprintf("%0.3f",auc(roc1))) )
}


legend("bottomright", aucText,lwd=2,bty="n",col=bioCol[1:(ncol(rt)-1)])

dev.off()

library(caret)
set.seed(123)
trainIndex <- createDataPartition(rt$Group, p = 0.5, list = FALSE)
test_data <- rt[trainIndex, ]
train_data <- rt[-trainIndex, ]
y=colnames(rt)[1]
bioCol=c("red","blue","green","yellow")
if(ncol(test_data)>4){
  bioCol=rainbow(ncol(test_data))}
pdf("test.pdf",width=5,height=5)
roc1=roc(test_data[,y], as.vector(test_data[,2]))
aucText=c( paste0(colnames(test_data)[2],", AUC=",sprintf("%0.3f",auc(roc1))) )
plot(roc1, col=bioCol[1])
for(i in 3:ncol(test_data)){
  roc1=roc(test_data[,y], as.vector(test_data[,i]))
  lines(roc1, col=bioCol[i-1])
  aucText=c(aucText, paste0(colnames(test_data)[i],", AUC=",sprintf("%0.3f",auc(roc1))) )
}
legend("bottomright", aucText,lwd=2,bty="n",col=bioCol[1:(ncol(test_data)-1)])
dev.off()

rt$Type <- factor(unlist(rt$Group), levels = c("n","o"))  # 或 rt$Type <- factor(rt$Type[,1])

vars <- c("TIMM50","RMND1","PNPO")
rt[vars] <- lapply(rt[vars], function(x) as.numeric(as.character(x)))

model <- glm(Type ~ TIMM50 + RMND1 + PNPO, data=rt, family="binomial")
rt$predicted <- predict(model, type="response")
roc_model <- roc(rt$Type, rt$predicted)
pdf("all.pdf",width=5,height=5)
plot(roc_model, col="red", lwd=2, main="Model ROC Curve")
text(0.5, 0.2, paste("AUC:", round(auc(roc_model), 3)), col="red")
ci_auc <- ci.auc(roc_model)
text(0.5, 0.15, paste("95% CI:", round(ci_auc[1], 3), "-", round(ci_auc[3], 3)), col="red")
dev.off()



test_data$Type <- factor(unlist(test_data$Group), levels = c("n","o")) 

vars <- c("TIMM50","RMND1","PNPO")
test_data[vars] <- lapply(test_data[vars], function(x) as.numeric(as.character(x)))

model <- glm(Type ~ TIMM50 + RMND1 + PNPO, data=test_data, family="binomial")
test_data$predicted <- predict(model, type="response")
roc_model <- roc(test_data$Type, test_data$predicted)
pdf("test.all.pdf",width=5,height=5)
plot(roc_model, col="red", lwd=2, main="Model ROC Curve")
text(0.5, 0.2, paste("AUC:", round(auc(roc_model), 3)), col="red")
ci_auc <- ci.auc(roc_model)
text(0.5, 0.15, paste("95% CI:", round(ci_auc[1], 3), "-", round(ci_auc[3], 3)), col="red")
dev.off()
