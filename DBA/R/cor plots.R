#' @import pheatmap
#' @import ape
cor.plots <- function(data.set, 
                      heatmap = TRUE, 
                      phylo = TRUE, 
                      method = "spearman", 
                      annotation = NULL, 
                      colors = NULL,
                      log = FALSE,
                      save = FALSE, wd = NULL){
  if (log == TRUE){
    cor.data <- cor(log(na.omit(data.set)+0.01), method = method)
  } else {
    cor.data <- cor(na.omit(data.set), method = method)
  }
  if (heatmap == TRUE){
    if (is.null(annotation) == TRUE){
      pheatmap(cor.data)
    } else {
      print(pheatmap(cor.data, annotation = annotation, annotation_colors = colors))
    }
  }
  if (phylo == TRUE){
    hc <- hclust(dist(1 - cor.data))
    plot.phylo(as.phylo(hc), type ="p", 
               edge.col = 4, edge.width = 3, 
               show.node.label = TRUE, 
               no.margin = TRUE)
  }
}

