

PCA <- function(dataset, PCA.Genes = FALSE){
  library(ggplot2)
  #generated PCA data
  genes.PCA <- prcomp(t(na.omit(expressed.genes)))
  
  #PCA variation 
  genes.PCA.var <- genes.PCA$sdev^2
  genes.PCA.var.per <- round(genes.PCA.var/sum(genes.PCA.var)*100, 1)
  
  #scree plots of variation 
  barplot(genes.PCA.var.per, main = "Scree Plot: Genes", 
          xlab = "Principal Component", 
          ylab = "Percent Variation")
  
  #organize data for ggplots
  genes.PCA.data <- data.frame(Sample = rownames(genes.PCA$x),
                               x = genes.PCA$x[,1],
                               y = genes.PCA$x[,2])
  
  #ggplot of PCA data
  print(ggplot(data = genes.PCA.data, aes(x = x, y = y, label = Sample))+
    geom_text(size = 2)+
    xlab(paste("PC1 - ", genes.PCA.var.per[1], "%", sep = ""))+
    ylab(paste("PC2 - ", genes.PCA.var.per[2], "%", sep = ""))+
    ggtitle("PCA: Expressed Genes"))
  
  # 100 genes the influence the PCA the greatest (either pos of neg)
  gene_score_ranked <- sort(abs(genes.PCA$rotation[,1]))
  gene_top_hun <- names(gene_score_ranked[1:100])
  genes.PCA$rotation[gene_top_hun, 1] # push to left on x axis (-) or right on x axis (+)
  
  if (PCA.Genes == TRUE) {
    return(gene_top_hun)
  }
  
}

genes <- PCA(expressed.genes, PCA.Genes = TRUE)
