#' @import pheatmap
#' @import ape
cor.plots <- function(data.set, 
                      heatmap = TRUE, 
                      phylo = TRUE, 
                      method = "spearman", 
                      annotation = NULL, 
                      colors = NULL,
                      save = FALSE, wd = NULL){
  
  cor.data <- cor(na.omit(data.set), method = method)
  if (save == TRUE){
    if(is.null(wd) == FALSE){
      setwd(wd)
    }
    par(mar=c(7,4,4,2)+0.1) 
    png(filename='Corr_heatmap.png', width=800, height=750)
  }
  if (heatmap == TRUE){
    if (is.null(annotation) == TRUE){
      pheatmap(cor.data)
    } else {
      pheatmap(cor.data, annotation = annotation, annotation_colors = colors)
    }
    if( save == TRUE){
      graphics.off()
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

