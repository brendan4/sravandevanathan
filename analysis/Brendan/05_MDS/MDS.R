
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)


# better plots *****NOTe IN THIS INSTANCE expressed.genes and expressed.trans 
#are without unmatched and are equiv to expressed.genes_WO22
par(mar=c(7,4,4,2)+0.1) 
png(filename='MDS Expressed Genes.png', width=750, height=600)
plotMDS(expressed.genes, cex = .75)
title("MDS Expresed Genes")
graphics.off()

par(mar=c(7,4,4,2)+0.1) 
png(filename='MDS Transcripts.png', width=750, height=600)
plotMDS(expressed.trans, cex = .75)
title("MDS Transcripts")
graphics.off()


