#' @import ggplot2
plot.var  <- function(dataset, 
                      zero.offset = 0.0001, 
                      title.name = "Gene Abundance"){
  
  dataset[is.na(dataset)] <- 0
  var <- apply(dataset,1, sd)
  plot(density(log(var + zero.offset)), 
       xlim=c(-15,10), 
       main = paste(title.name,
                    ": Density x = log(var + zero.offset)", 
                    sep = ""))
}






