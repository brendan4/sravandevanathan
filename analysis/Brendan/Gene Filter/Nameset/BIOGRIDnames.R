library(dplyr)
gene.list <- read.delim("BIOGRIDnames.txt", sep = "\t")
gene.list <- gene.list[,c(8,13)]
gene.list <- unique.data.frame(gene.list)
gene.list <- t(gene.list)
gene.names <- gene.list[1,]
filt.data <- filter.genes(expressed.genes, gene.names)
# hierarchical clustering of samples and heatmap of sample similarities

ann_colors = list(Var1 = Var1)
library(pheatmap)
pheatmap(cor(na.omit(filt.data), method = "spearman")) 

#simple clustering based on FPKMs
d <- cor(na.omit(filt.data), method="spearman")
hc <- hclust(dist(1-d))
plot.phylo(as.phylo(hc), type="p", edge.col=4, edge.width=3, show.node.label=TRUE, no.margin=TRUE)
dev.off()