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


# better plots *****NOTe IN THIS INSTANCE expressed.genes and expressed.trans 
#are without unmatched and are equiv to expressed.genes_WO22
par(mar=c(7,4,4,2)+0.1) 
png(filename='MDS Expressed Genes - improved.png', width=750, height=600)
plotMDS(expressed.genes, cex = .75)
title("MDS Expresed Genes")
graphics.off()

par(mar=c(7,4,4,2)+0.1) 
png(filename='MDS Trans Data - improved.png', width=750, height=600)
plotMDS(expressed.trans, cex = .75)
title("MDS t_data")
graphics.off()


