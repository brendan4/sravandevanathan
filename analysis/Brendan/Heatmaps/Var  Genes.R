library(gplots)
library(ComplexHeatmap)

# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(log2(na.omit(expressed.genes)), 1, var)
head(var_genes)

# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:100]
head(select_var)

# Subset logcounts matrix
highly_variable <- expressed.genes[select_var,]
dim(highly_variable)


# We can split the heatmap into clusters
fontsize <- 0.6

Heatmap(as.matrix(log10(highly_variable)),
        cluster_columns=FALSE,
        row_names_side = "left",
        row_hclust_side = "left",
        row_names_gp=gpar(cex=fontsize),
        row_hclust_width = unit(3, "cm"),
        clustering_distance_rows ="maximum",
        clustering_method_rows = "centroid",
        km=2) # number of clusters you want

library(ggplot2)

#princple component for expressed.genes
genes.PCA <- prcomp(t(na.omit(highly_variable)), scale. = T)
summary(genes.PCA)

#precent variance 
genes.PCA.var <- genes.PCA$sdev^2
genes.PCA.var.per <- round(genes.PCA.var/sum(genes.PCA.var)*100, 1)

#scree plots of variation 
barplot(genes.PCA.var.per, main = "Scree Plot: Genes", 
        xlab = "Principal Component", 
        ylab = "Percent Variation")

# reformate of PCA data into a data frame for ggplot usage
genes.PCA.data <- data.frame(Sample = rownames(genes.PCA$x),
                             x = genes.PCA$x[,1],
                             y = genes.PCA$x[,2])

#ggplot of PCA data
ggplot(data = genes.PCA.data, aes(x = x, y = y, label = Sample))+
  geom_text(size = 2)+
  xlab(paste("PC1 - ", genes.PCA.var.per[1], "%", sep = ""))+
  ylab(paste("PC2 - ", genes.PCA.var.per[2], "%", sep = ""))+
  ggtitle("PCA: Expressed Genes")

#gene_PCA_plot <- ggbiplot(genes.PCA)

# 100 genes the influence the PCA the greatest (either pos of neg)
gene_score_ranked <- sort(abs(genes.PCA$rotation[,1]))
gene_top_hun <- names(gene_score_ranked[1:100])
genes.PCA$rotation[gene_top_hun, 1] # push to left on x axis (-) or right on x axis (+)

