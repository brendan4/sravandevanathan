#' @import pheatmap
#' @import ape
cor.plots <- function(data.set, heatmap = TRUE, phylo = TRUE, method = "spearman"){
  
  cor.data <- cor(na.omit(data.set), method = method)
  if (heatmap == TRUE){
    pheatmap(cor.data)
  }
  if (phylo == TRUE){
    hc <- hclust(dist(1 - cor.data))
    plot.phylo(as.phylo(hc), type ="p", 
               edge.col = 4, edge.width = 3, 
               show.node.label = TRUE, 
               no.margin = TRUE)
  }
}