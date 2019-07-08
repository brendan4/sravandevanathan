#Uses data generated from PCA.R scripts: gene_top_hun and trans_top_hun

genes_top <- filter.genes(expressed.genes, gene_top_hun)
trans_top <- filter.genes(expressed.trans, trans_top_hun)

boxplot(log2(genes_top+0.1), 
        names=colnames(genes_top), las=2, ylab="log2(FPKM)", 
        main="Distribution of FPKMs for all libraries")
boxplot(log2(trans_top+0.1), 
        names=colnames(trans_top), las=2, ylab="log2(FPKM)", 
        main="Distribution of FPKMs for all libraries")

# color and packages 
require("RColorBrewer")
library(zFPKM)
require("gplots")
myCol <- colorRampPalette(c("dodgerblue", "black", "yellow"))(100)
myBreaks <- seq(-3, 3, length.out=101)

heat <- zFPKM(trans_top)  ####swtiched here for trans vs genes
heat <- heat[is.finite(rowSums(heat)),] # filters out -inf values 

par(mar=c(7,4,4,2)+0.1) 
png(filename='Euclidean PCA geneset.png', width=800, height=750)

#Euclidean distance
heatmap.2(as.matrix(heat),
          col=myCol,
          breaks=myBreaks,
          main="Euclidean Distance",
          key=T, keysize=.5,
          scale="none",
          density.info="none",
          reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
          trace="none",
          cexRow=.2,
          cexCol=0.8,
          distfun=function(x) dist(x, method="euclidean"),
          hclustfun=function(x) hclust(x, method="average"))
graphics.off()


par(mar=c(7,4,4,2)+0.1) 
png(filename='Pearson correlation PCA geneset.png', width=800, height=750)

#1-cor distance
heatmap.2(as.matrix(heat),
          col=myCol,
          breaks=myBreaks,
          main="Person Correlation",
          key=T, keysize=1.0,
          scale="none",
          density.info="none",
          reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
          trace="none",
          cexRow=0.2,
          cexCol=0.8,
          distfun=function(x) as.dist(1-cor(t(x))),
          hclustfun=function(x) hclust(x, method="average"))
graphics.off()

rm(heat,logcounts_genes,myBreaks,myCol, var_genes)
rm(select_var_genes)
