#' @import ggplot2
MDS <- function(data.set, scaled = FALSE){
  
  #distance - Euclidean 
  if (scaled == TRUE){
    d_genes <- dist(scale(t(data.set),
                          center = TRUE, 
                          scale = TRUE), 
                    method = "euclidean") # for scaled values
  }else {
    d_genes <- dist(t(data.set), method = "euclidean")
  }
  
  #MDS with euclidean distance
  mds_genes <- cmdscale(d_genes, eig = TRUE, x.ret = TRUE)
  
  #variation percentage 
  mds.var.per_genes <- round(mds_genes$eig/sum(mds_genes$eig)*100,1)
  
  #organization of data for ggplot use
  mds.values_genes  <- mds_genes$points
  mds.data_genes <- data.frame(Sample = rownames(mds.values_genes),
                               X = mds.values_genes[,1],
                               Y = mds.values_genes[,2])
  
  #plot of MDS
  ggplot(data = mds.data_genes, aes(x=X, y=Y, label = Sample))+
    geom_text(size = 2.5)+
    theme_bw()+
    xlab(paste("MDS1 - ", mds.var.per_genes[1], "%", sep = ""))+
    ylab(paste("MD2 - ", mds.var.per_genes[2], "%", sep = ""))+
    ggtitle("MDS Euclidean distance - Genes")
  
}
  