#princple component for filtered data: expressed.genes and expressed.trans
rownames(expressed.genes) <- expressed.genes[,1]
rownames(expressed.trans) <- expressed.trans[,1]
PCA_expressed.genes <- expressed.genes[-1]
PCA_expressed.trans <- expressed.trans[-1]
NA_genes<- PCA_expressed.genes[is.na(PCA_expressed.genes),]

library(ggplot2)
library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)


#princple component for expressed.genes
genes.PCA <- prcomp(t(log(na.omit(PCA_expressed.genes) +0.01)),scale = TRUE)
summary(genes.PCA)

#basic plot of PCA1 and PCA 2
plot(genes.PCA$x[,1], genes.PCA$x[,2],
     main = "PCA: Genes - log scale",
     xlab = paste("PC1 - ", genes.PCA.var.per[1], "%", sep = ""), 
     ylab = paste("PC2 - ", genes.PCA.var.per[2], "%", sep = ""))

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



#princple component for t_data
trans.PCA <- prcomp(t(na.omit(PCA_expressed.trans)),scale = TRUE)
summary(trans.PCA)

plot(trans.PCA$x[,1], trans.PCA$x[,2],
     main = "PCA: t_data",
     xlab = paste("PC1 - ", trans.PCA.var.per[1], "%", sep = ""), 
     ylab = paste("PC2 - ", trans.PCA.var.per[2], "%", sep = ""))

trans.PCA.var <- trans.PCA$sdev^2
trans.PCA.var.per <- round(trans.PCA.var/sum(trans.PCA.var)*100, 1)

barplot(trans.PCA.var.per, main = "Scree Plot: Trans", 
        xlab = "Principal Component", 
        ylab = "Percent Variation")

trans.PCA.data <- data.frame(Sample = rownames(trans.PCA$x),
                             x = trans.PCA$x[,1],
                             y = trans.PCA$x[,2])

ggplot(data = trans.PCA.data, aes(x = x, y = y, label = Sample))+
  geom_text(size = 1.5)+
  xlab(paste("PC1 - ", trans.PCA.var.per[1], "%", sep = ""))+
  ylab(paste("PC2 - ", trans.PCA.var.per[2], "%", sep = ""))+
  ggtitle("PCA: t_data")

trans_score_ranked <- sort(abs(trans.PCA$rotation[,1]))
trans_top_hun <- names(trans_score_ranked[1:100])
trans.PCA$rotation[trans_top_hun, 1] # push to left on x axis (-) or right on x axis (+)

#t( log(expressed.genes[-1]  + 0.01)  )

