library(randomForest)
set.seed(123)
inputFile <- "diffGeneExp.txt"
data <- read.table(inputFile, header = TRUE, sep = "\t",
                   check.names = FALSE, row.names = 1)
data <- t(data)
group <- gsub("(.*)\\_(.*)", "\\2", rownames(data))
group <- factor(group)
vars <- apply(data, 2, function(x) var(as.numeric(x)) > 0)
data <- data[, vars, drop = FALSE]
cls_tab <- table(group)
min_n   <- min(cls_tab)
sampv   <- rep(min_n, length(cls_tab))  # 按 group 的水平顺序
p <- ncol(data)
set.seed(123)
rf <- randomForest(as.factor(group) ~ ., data = data,
                   ntree = 1500,
                   mtry  = max(2, floor(p/3)),  # 可按需调大或用 tuneRF()
                   sampsize = sampv,
                   importance = TRUE)
pdf("forest.pdf", width = 6, height = 6)
plot(rf, main = "Random forest (balanced, mtry=p/3, ntree=1500)", lwd = 2)
dev.off()
imp_gini <- importance(rf, type = 2, scale = FALSE)  # 非负，很多为0是正常
write.table(data.frame(Gene = rownames(imp_gini),
                       MeanDecreaseGini = imp_gini[, "MeanDecreaseGini"]),
                       "importance_gini.txt", sep = "\t", quote = FALSE, row.names = FALSE)
TopN <- 10
imp_df <- data.frame(Gene = rownames(imp_gini),
                     Score = imp_gini[, "MeanDecreaseGini"])
imp_df <- imp_df[order(imp_df$Score, decreasing = TRUE), ]
top_df <- head(imp_df, TopN)
suppressPackageStartupMessages({library(ggplot2); library(scales) })
pdf("geneImportance_topN.pdf", width = 7, height = 6)
ggplot(top_df, aes(x = Score, y = reorder(Gene, Score))) +
               geom_segment(aes(x = 0, xend = Score, yend = reorder(Gene, Score)), linewidth = 0.6) +
               geom_point(size = 1.8) +
               labs(x = "MeanDecreaseGini (Top N)", y = NULL) +
               scale_x_continuous(breaks = pretty_breaks(8), expand = expansion(mult = c(0, 0.02))) +
               theme_bw(base_size = 11) +
               theme(panel.grid.minor = element_blank(),
               axis.text.y = element_text(size = 8))
 dev.off()

rfGenes=read.table("importance_gini.txt", header = TRUE, sep = "\t",
                   check.names = FALSE, row.names = 1)
top5_names <- rownames(rfGenes)[order(rfGenes$MeanDecreaseGini, decreasing = TRUE, na.last = NA)][1:5]
write.table(top5_names, file="rfGenes.txt", sep="\t", quote=F, col.names=F, row.names=F)

sigExp=t(data[,top5_names])
sigExpOut=rbind(ID=colnames(sigExp),sigExp)
write.table(sigExpOut, file="rfGeneExp.txt", sep="\t", quote=F, col.names=F)


