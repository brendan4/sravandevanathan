library(dplyr)

pData(bg_filt) = data.frame(id = sampleNames(bg_filt), group = pheno)
gene.results = stattest(bg_filt, feature="gene", covariate="group.pheno", getFC=TRUE, meas="FPKM")

#gene names 
indices <- match(gene.results$id, texpr(bg_filt, 'all')$gene_id)
gene_names_for_result <- texpr(bg_filt, 'all')$gene_name[indices]
gene.results <- data.frame(geneNames = gene_names_for_result, gene.results)

#processing 
gene.results = arrange(gene.results, pval)
filt.names <- subset(gene.results, gene.results$qval<0.05)

#filter out . from data set as it acts as a regular expression 
filt.set <- filter.out.genes(filt.names, gene.list = c("\\."), by.rownames = FALSE, col = 1)
diff.genes <- filter.genes(expressed.genes, filt.set$geneNames)

genesPCA1 <- PCA(diff.genes, PCA.Genes = TRUE, pheno = pheno)

PCA.filt <- filter.genes(diff.genes, gene.list = genesPCA1)

library(gplots)
library(ComplexHeatmap)

PCA.filt <- pretty.gene.name(PCA.filt, as.row.names = TRUE, remove.dups = TRUE)

Heatmap(t(as.matrix(log2(PCA.filt+0.01))),
        column_names_side = "bottom",
        row_names_side = "left",
        row_hclust_side = "left",
        row_names_gp=gpar(cex=0.6),
        row_hclust_width = unit(3, "cm"),
        clustering_distance_rows ="euclidean",
        clustering_method_rows = "centroid",
        km=3) # number of clusters you want
graphics.off()
