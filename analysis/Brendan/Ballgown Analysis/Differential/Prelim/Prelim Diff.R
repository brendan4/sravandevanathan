library(dplyr)
data(pheno)

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

#PCA with diff genes
genesPCA1 <- PCA(diff.genes, PCA.Genes = TRUE, pheno = pheno)

# pca one genes 
PCA.filt <- filter.genes(diff.genes, gene.list = genesPCA1)

library(gplots)
library(ComplexHeatmap)

PCA.filt <- pretty.gene.name(PCA.filt, as.row.names = TRUE, remove.dups = TRUE)

#heatmaps with PCA1 genes
pheno.types <- t(pheno[,-3])
pheno.types <- t(pheno)
pheno.types <- pheno.types[-1,]
names(pheno.types) <- t(pheno[,1])

#saving options 
par(mar=c(7,4,4,2)+0.1) 
png(filename='Heatmap - PCA 1.png', width=800, height=750)

Heatmap(as.matrix(log2(PCA.filt+0.01)),
        column_names_side = "top",
        row_names_side = "left",
        row_hclust_side = "left",
        row_names_gp=gpar(cex=0.6),
        row_hclust_width = unit(3, "cm"),
        clustering_distance_rows ="euclidean",
        clustering_method_rows = "centroid",
        top_annotation = HeatmapAnnotation(as.data.frame(pheno.types),
                                              which = "column",
                                              show_legend=TRUE))
graphics.off()

