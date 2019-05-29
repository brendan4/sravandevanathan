# color and packages 
require("RColorBrewer")
library(zFPKM)
require("gplots")
myCol <- colorRampPalette(c("dodgerblue", "black", "yellow"))(100)
myBreaks <- seq(-3, 3, length.out=101)

# Get log2 counts per million
logcounts_genes <- cpm(na.omit(expressed.genes_WO22),log=TRUE)
# estimate the variance for each row in the logcounts matrix
var_genes <- apply(logcounts_genes, 1, var)
# Get the gene names for the top 100 most variable genes
select_var_genes <- names(sort(var_genes, decreasing=TRUE))[1:100]
heat <- zFPKM(na.omit(expressed.genes_WO22))  
heat <- heat[which(rownames(heat) %in% select_var_genes), ] # using highly expressed genes to filter heat
heat <- heat[is.finite(rowSums(heat)),] # filters out -inf values 


par(mar=c(7,4,4,2)+0.1) 
png(filename='Euclidean_Distance.png', width=800, height=750)

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
          hclustfun=function(x) hclust(x, method="ward.D2"))
graphics.off()


par(mar=c(7,4,4,2)+0.1) 
png(filename='Pearson correlation.png', width=800, height=750)

#1-cor distance
heatmap.2(as.matrix(heat),
          col=myCol,
          breaks=myBreaks,
          main="Title",
          key=T, keysize=1.0,
          scale="none",
          density.info="none",
          reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
          trace="none",
          cexRow=0.2,
          cexCol=0.8,
          distfun=function(x) as.dist(1-cor(t(x))),
          hclustfun=function(x) hclust(x, method="ward.D2"))
graphics.off()

