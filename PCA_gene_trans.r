#princple component for filtered data: expressed.genes and expressed.trans
rownames(expressed.genes) <- expressed.genes[,1]
rownames(expressed.trans) <- expressed.trans[,1]
PCA_expressed.genes <- expressed.genes[-1]
PCA_expressed.trans <- expressed.trans[-1]
NA_genes<- PCA_expressed.genes[is.na(PCA_expressed.genes),]

library(ggplot2)

#princple component for expressed.genes
genes.PCA <- prcomp(t(na.omit(PCA_expressed.genes)),scale = TRUE)
summary(genes.PCA)
plot(genes.PCA$x[,1], genes.PCA$x[,2])
genes.PCA.var <- genes.PCA$sdev^2
genes.PCA.var.per <- round(genes.PCA.var/sum(genes.PCA.var)*100, 1)
barplot(genes.PCA.var.per, main = "Scree Plot: Gene", xlab = "Principal Component", ylab = "Percent Variation")
genes.PCA.data <- data.frame(Sample = rownames(genes.PCA$x),
                             x = genes.PCA$x[,1],
                             y = genes.PCA$x[,2])
ggplot(data = genes.PCA.data, aes(x = x, y = y, label = Sample))+
  geom_text(size = 3)+
  xlab(paste("PC1 - ", genes.PCA.var.per[1], "%", sep = ""))+
  ylab(paste("PC2 - ", genes.PCA.var.per[2], "%", sep = ""))+
  ggtitle("PCA: Expressed Genes")
gene_score_ranked <- sort(abs(genes.PCA$rotation[,1]))
gene_top_ten <- names(gene_score_ranked[1:10])
genes.PCA$rotation[gene_top_ten, 1] # push to left on x axis (-) or right on x axis (+)



#princple component for t_data
trans.PCA <- prcomp(t(na.omit(PCA_expressed.trans)),scale = TRUE)
summary(trans.PCA)
plot(trans.PCA$x[,1], trans.PCA$x[,2])
trans.PCA.var <- trans.PCA$sdev^2
trans.PCA.var.per <- round(trans.PCA.var/sum(trans.PCA.var)*100, 1)
barplot(trans.PCA.var.per, main = "Scree Plot: Trans", xlab = "Principal Component", ylab = "Percent Variation")
trans.PCA.data <- data.frame(Sample = rownames(trans.PCA$x),
                             x = trans.PCA$x[,1],
                             y = trans.PCA$x[,2])
ggplot(data = trans.PCA.data, aes(x = x, y = y, label = Sample))+
  geom_text(size = 2)+
  xlab(paste("PC1 - ", trans.PCA.var.per[1], "%", sep = ""))+
  ylab(paste("PC2 - ", trans.PCA.var.per[2], "%", sep = ""))+
  ggtitle("PCA: t_data")
trans_score_ranked <- sort(abs(trans.PCA$rotation[,1]))
trans_top_ten <- names(trans_score_ranked[1:10])
trans.PCA$rotation[trans_top_ten, 1] # push to left on x axis (-) or right on x axis (+)
