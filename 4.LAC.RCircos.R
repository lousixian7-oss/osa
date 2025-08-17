library(tidyverse)
REFFile = "REF.txt"
GeneFile = "diffgene.txt"
REF=read.table(REFFile, header=T, sep="\t")
Gene=read.table(GeneFile, header=T, sep="\t")
a=merge(Gene,REF,by="Gene")
a$Chr <- paste('chr', a$Chr,sep = "")
library("RCircos")
cytoBandIdeogram=read.table("refer.txt", header=T, sep="\t", check.names=F)
chr.exclude <- NULL
cyto.info <- cytoBandIdeogram
tracks.inside <- 4
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)
rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.params$text.size=0.8
rcircos.params$point.size=5
RCircos.Reset.Plot.Parameters(rcircos.params)
pdf(file="RCircos.pdf", width=11, height=11)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
name.col <- 4
side <- "in"
track.num <- 1
a <- a[, c(2, 3, 4, 1, if (ncol(a) > 4) 5:ncol(a))]
RCircos.Gene.Connector.Plot(a, track.num, side)
track.num <- 2
RCircos.Gene.Name.Plot(a, name.col, track.num, side)
dev.off()
