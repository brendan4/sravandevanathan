# Get log2 counts per million
logcounts_genes <- cpm(na.omit(expressed.genes_WO22),log=TRUE)
boxplot(logcounts_genes, xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(logcounts_genes),col="blue")
title("Boxplots of logCPMs (unnormalised)")

# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(logcounts_genes, 1, var)
# Get the gene names for the top 500 most variable genes
select_var_genes <- names(sort(var_genes, decreasing=TRUE))[1:500]
# Subset logcounts matrix
highly_variable_lcpm_genes <- logcounts_genes[select_var_genes,]
dim(highly_variable_lcpm_genes)

mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
heatmap.2(highly_variable_lcpm_genes,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",scale="row")
