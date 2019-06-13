library(gplots)
library(ComplexHeatmap)


#estimate the variance for each row in the logcounts matrix
var_genes <- apply(log2(na.omit(expressed.genes)+0.1), 1, var)
head(var_genes)

# Get the gene names for the top 100 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:100]
head(select_var)

# Subset logcounts matrix
highly_variable <- expressed.genes[select_var,]
dim(highly_variable)

highly_variable$var<- apply(log2(highly_variable+0.1), 1, var)
head(var_genes)

###################SKIP: if don't want sex related genes filtered out######################
#sex insight
sex <- t(sex)
sex <- highly_variable[grep("^XIST", rownames(highly_variable)),]
sex[,which(sex > 10)] = "F"
sex[,which(sex < 10)] = "M"
rownames(sex) <- "Gender"
sex <- t(sex)

#bg_filt diff expression with sex as covariate: uses fitler ot func in gene filter/ tools

sex.filt <- filter.out.genes(expressed.genes, filt.names$geneNames)
sex.filt <- na.omit(sex.filt)


# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(log2(sex.filt), 1, var)
head(var_genes)

# Get the gene names for the top 100 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:100]
head(select_var)

# Subset logcounts matrix
highly_variable <- expressed.genes[select_var,]
dim(highly_variable)

#OPTIONAL: make gene names pretty: load pretty func in brendan/tools 
highly_variable <- pretty.gene.name(highly_variable) # make the pretty col 
highly_variable <- highly_variable[-which(duplicated(highly_variable$pretty)),] # drop the duplicated 
rownames(highly_variable) <- highly_variable$pretty # make gene names pretty
highly_variable <- highly_variable[,-which(colnames(highly_variable) %in% c("pretty"))] # drop pretty col

fontsize <- 0.6

par(mar=c(7,4,4,2)+0.1) 
png(filename='Sex_filtered.png', width=800, height=750)

# use t() for k means clustering of sample no genes 
Heatmap(t(as.matrix(log2(highly_variable))),
        column_names_side = "bottom",
        row_names_side = "left",
        row_hclust_side = "left",
        row_names_gp=gpar(cex=fontsize),
        row_hclust_width = unit(3, "cm"),
        clustering_distance_rows ="manhattan",
        clustering_method_rows = "centroid",
        km=3) # number of clusters you want
graphics.off()


#simple PCA plots

library(ggplot2)

#princple component for expressed.genes
genes.PCA <- prcomp(t(na.omit(highly_variable)))
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
gene_top_hun <- names(gene_score_ranked[1:10])
genes.PCA$rotation[gene_top_hun, 1] # push to left on x axis (-) or right on x axis (+)

