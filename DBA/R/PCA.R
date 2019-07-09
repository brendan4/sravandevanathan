PCA <- function(dataset, scaled = FALSE, 
                PCA.Genes = FALSE, 
                pheno = NULL, 
                label.size = 2, 
                pca.dim = c(1,2)){
  
  #generated PCA data
  if (scaled == TRUE){
    genes.PCA <- prcomp(log(t(na.omit(dataset))+0.01), scale. = TRUE)
  }else {
    genes.PCA <- prcomp(t(na.omit(dataset)))
  }
  
  #PCA variation 
  genes.PCA.var <- genes.PCA$sdev^2
  genes.PCA.var.per <- round(genes.PCA.var/sum(genes.PCA.var)*100, 1)
  
  #scree plots of variation 
  barplot(genes.PCA.var.per, main = "Scree Plot: Genes", 
          xlab = "Principal Component", 
          ylab = "Percent Variation")
  
  #organize data for ggplots
  genes.PCA.data <- data.frame(Sample = rownames(genes.PCA$x),
                               x = genes.PCA$x[,pca.dim[1]],
                               y = genes.PCA$x[,pca.dim[2]])
  
  if (is.null(pheno) == TRUE){
  #ggplot of PCA data
  print(ggplot(data = genes.PCA.data, aes(x = x, y = y, label = Sample))+
          geom_text(size = label.size)+
          xlab(paste("PC1 - ", genes.PCA.var.per[pca.dim[1]], "%", sep = ""))+
          ylab(paste("PC2 - ", genes.PCA.var.per[pca.dim[2]], "%", sep = ""))+
          ggtitle("PCA: Expressed Genes"))
  }else {
    #ggplot of PCA data
    colors <- c("green", "orange", "red", "blue")
    names(colors) <- c(pheno$pheno[!duplicated(pheno$pheno)])
    print(ggplot(data = genes.PCA.data, aes(x = x, y = y, label = Sample))+
            geom_text(size = label.size, aes(color = pheno$pheno))+
            scale_color_manual(breaks = c("8", "6", "4", "2"),
                               values=colors) +
            xlab(paste("PC1 - ", genes.PCA.var.per[pca.dim[1]], "%", sep = ""))+
            ylab(paste("PC2 - ", genes.PCA.var.per[pca.dim[2]], "%", sep = ""))+
            ggtitle("PCA: Expressed Genes"))
  }
  
  # 100 genes the influence the PCA the greatest (either pos of neg)
  gene_score_ranked <- sort(abs(genes.PCA$rotation[,pca.dim[1]]))
  gene_top_hun <- names(gene_score_ranked[1:100])
  genes.PCA$rotation[gene_top_hun, 1] # push to left on x axis (-) or right on x axis (+)
  
  if (PCA.Genes == TRUE) {
    return(gene_top_hun)
  }
  
}