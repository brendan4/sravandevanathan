#' Gene Variation Between Rows
#' 
#' Filters genes non-specficically from a data table and calculates variation between samples. Returns a dataset
#' with a var col.

#' @param data.set Data.frame with rownames as genes and cols as different samples.
#' @param gene.list A vector with gene names or name to be filtred 
#' @param graph a bolean stating if a simple bar graph should be returned 
#' @param pretty.names a bolean stating if gene names should returned without unique start and stop identifier:
#' may create duplicates. In such a case a pretty col will be returned instead. See pretty.gene.names func.
#' @param func which log scale should the data be  subjected to
#' @param zero addition to offset taking log of zero. Otherwise zeros in the dataset will create inf that mutate entire roww
#' 
#' @return Returns a new data frame with only genes matching to genes.list and a variation column.
#' If a gene.list name fails to match it will be stated. See gene.filter func. 
#' If pretty.names is TRUE and duplicates exist a pretty gene col will be return in the data.frame
#' 
#' 
#' @examples 
#' filtered.data <- Var.samples(expressed.genes, gene.list = c("RPL11"))
#' filtered.data <- Var.samples(expressed.genes, gene.list = c("RPL11"), pretty.names = TRUE, graph = TRUE)

var.samples <- function(data.set, gene.list ,graph = FALSE, pretty.names = FALSE, func = log2, zero = 0.01){
  
  filt.data <- filter.genes(data.set, gene.list = gene.list)
  filt.data$var <- apply(func(filt.data + zero),1, var)
  
  if (pretty.names == TRUE){
     filt.data <- pretty.gene.name(filt.data)
     if (sum(duplicated(filt.data$Pretty)) == 0){
       rownames(filt.data) <- filt.data$pretty
       filt.data <- filt.data[,-which(colnames(filt.data) %in% c("pretty"))]
     }else{
       print("Warning: Duplicates Found")
     }
  }
  
  if (graph == TRUE){
    print(ggplot(filt.data, aes(x = rownames(filt.data) , y = var)) + geom_bar(stat="identity", fill="tomato3"))
  }
  return(filt.data)
}
