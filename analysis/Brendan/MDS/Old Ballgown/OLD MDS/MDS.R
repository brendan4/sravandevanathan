#setting rows names and removing names col 
rownames(expressed.genes) <- expressed.genes[,1]
rownames(expressed.trans) <- expressed.trans[,1]
PCA_expressed.genes <- expressed.genes[-1]
PCA_expressed.trans <- expressed.trans[-1]

library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)

#MDS on expressed genes with unmatched 
plotMDS(PCA_expressed.genes)
title("MDS Expressed Genes")

#MDS on t_data with unmatched 
plotMDS(PCA_expressed.trans)
title("MDS t_data")

#removed unmatched 
expressed.genes_WO22 <- PCA_expressed.genes[-22]
plotMDS(expressed.genes_WO22)
title("MDS Expresed Genes WO unmatched")

expressed.trans_WO22 <- PCA_expressed.trans[-22]
plotMDS(expressed.trans_WO22)
title("MDS t_data WO unmatched")
