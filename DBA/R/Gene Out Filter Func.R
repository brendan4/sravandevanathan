#' Filter out genes 
#' 
#' Filters out genes non-specficically from a data table

#' @param data.table Data.frame with rownames as genes.
#' @param gene.list A vector with gene names or name to be filtred out
#' @param by.rownames default filters by row names of the data.table. If false must provide col index.
#' @param col Index of column postion if rownames are not to be used. Must provide col index if by.row is false 
#' 
#' @return Returns a new data frame with genes matching to genes.list filtered out.
#' 
#' @examples 
#' filter.out.genes(expressed.genes, c("RPL11", "RPS", by.rownames = TRUE))
#' filter.out.genes(expressed.genes, c("RP"), byrownames = FALSE, col = 2)


#filers data set with genelist using rownames 
filter.out.genes <- function(data.table, gene.list, by.rownames = TRUE, col = NULL){
  filtered.data <- data.table
  for (name in gene.list){
    print(name)
    
    if(by.rownames == TRUE){
      info <- grepl(paste("^",name, sep = ""),
                    rownames(filtered.data))
    }else {
      if (is.null(col) == TRUE){
        stop("No column index given and by.rownames is FALSE")
      }else{
      info <- grepl(paste("^",name, sep = ""),
                    filtered.data[ ,col])
      }
    }
    filtered.data <- filtered.data[which(info == F),] #matches pattern to dataframe and index
  }
  return(filtered.data)
}

