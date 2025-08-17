
library(ConsensusClusterPlus)   
expFile="diffGeneExp.txt"         
workDir="C:\\Users\\23634\\Desktop\\论文代码"

data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.matrix(data)

group=sapply(strsplit(colnames(data),"\\_"), "[", 2)
data=data[,group=="o"]

storage.mode(data) <- "numeric"
gene_sd <- apply(data, 1, sd, na.rm = TRUE)
data <- data[gene_sd > 1e-8, , drop = FALSE]

maxK=9
results=ConsensusClusterPlus(data,
              maxK=maxK,
              reps=500,
              pItem=0.8,
              pFeature=1,
              title=workDir,
              clusterAlg="pam",
              distance="euclidean",
              seed=123456,
              plot="png")


clusterNum=3       
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("LACcluster")
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(cluster$LACcluster))
cluster$LACcluster=letter[match(cluster$LACcluster, uniqClu)]
outTab=cbind(t(data), cluster)
outTab=rbind(ID=colnames(outTab), outTab)
write.table(outTab, file="LACcluster.txt", sep="\t", quote=F, col.names=F)


  library(pheatmap)
  library(cluster)
  library(GSVA)    


## ===== 1) Consensus matrix heatmap + consensus scores =====
k_use <- clusterNum
res_k <- results[[k_use]]

# Get consensus matrix (different CCP versions use different slots)
cm <- if (!is.null(res_k$consensusMatrix)) res_k$consensusMatrix else res_k$ml
stopifnot(is.matrix(cm), nrow(cm) > 1)

# Cluster labels
cls <- res_k$consensusClass

# Ensure dimnames exist and align with class names
if (is.null(colnames(cm))) {
  nm <- names(cls); if (is.null(nm)) nm <- paste0("S", seq_len(ncol(cm)))
  rownames(cm) <- colnames(cm) <- nm
}
if (is.null(names(cls))) names(cls) <- colnames(cm)

# Sample-wise consensus score (mean within-cluster consensus excluding self)
intra_cons <- sapply(names(cls), function(s){
  ix <- names(cls)[cls == cls[s]]
  v  <- cm[s, ix]
  if (length(v) > 1) mean(v[v < 0.999999], na.rm = TRUE) else NA_real_
})

# Cluster-level mean consensus (upper triangle mean)
cluster_cons <- tapply(names(cls), cls, function(ix){
  m  <- cm[ix, ix, drop = FALSE]
  ut <- m[upper.tri(m)]
  if (length(ut)) mean(ut, na.rm = TRUE) else NA_real_
})

# PAC (proportion of ambiguous clustering; lower is better)
u <- cm[upper.tri(cm)]
u <- u[is.finite(u) & !is.na(u)]
CDFx <- function(x, v) mean(v <= x, na.rm = TRUE)
PAC  <- CDFx(0.9, u) - CDFx(0.1, u)

# Mean silhouette (higher is better), distance = 1 - consensus
d <- as.dist(1 - cm)
sil <- silhouette(as.integer(factor(cls)), d)
mean_sil <- mean(sil[, "sil_width"])

# Annotation for heatmap
ann <- data.frame(
  Cluster = factor(paste0("C", cls), levels = paste0("C", sort(unique(cls)))),
  ConsensusScore = intra_cons[colnames(cm)]
)
rownames(ann) <- colnames(cm)
ann_colors <- list(
  Cluster = setNames(RColorBrewer::brewer.pal(n = length(levels(ann$Cluster)), "Set1"),
                     levels(ann$Cluster)),
  ConsensusScore = colorRampPalette(c("#f7fbff", "#08306b"))(100)
)

pdf("consensus_matrix_k_selected.pdf", width = 6.5, height = 6.2)
pheatmap(cm,
         annotation_col = ann,
         annotation_colors = ann_colors,
         color = colorRampPalette(c("white", "royalblue"))(100),
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = FALSE, show_colnames = FALSE,
         main = sprintf("Consensus matrix (k=%d)\nPAC=%.3f, mean silhouette=%.3f",
                        k_use, PAC, mean_sil)
)
dev.off()

write.table(data.frame(
  k = k_use, PAC = round(PAC,3), mean_silhouette = round(mean_sil,3),
  t(round(unlist(cluster_cons),3))),
  "consensus_stability_metrics.txt", sep = "\t", quote = FALSE, row.names = FALSE)

