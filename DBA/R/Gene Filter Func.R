#' Filter genes 
#' 
#' Filters genes non-specficically from a data table and returns that subset

#' @param data.table Data.frame with gene names
#' @param gene.list A vector with gene names or name to be filtred 
#' @param by.rownames default filters by row names of the data.table. If false must provide col index.
#' @param col Index of column postion if rownames are not to be used. Must provide col index if by.row is false 
#' @param lazy Should function filter by pretty gene names. Works only with rownames
#' 
#' @return Returns a new data frame with only genes matching to genes.list. 
#' If a gene.list name fails to match it will be stated.
#' 
#' @examples 
#' filter.genes(expressed.genes, c("RPL11", "RPS", by.rownames = TRUE))
#' filter.genes(expressed.genes, c("RP"), byrownames = FALSE, col = 2)
#' diff.genes <- filter.genes(na.omit(expressed.genes), filt.names$geneNames)


#filers data set with genelist using rownames 
filter.genes <- function(data.table, gene.list, by.rownames = TRUE, col = NULL, lazy = TRUE){
  filtered.data <- data.frame()
  if (lazy == FALSE){
    for(row in 1:nrow(data.table)){
      new.name <- str_extract(rownames(data.table[row,]), "^[-a-zA-Z0-9]*")
      data.table$pretty[row] <- new.name 
    }
  }
  for (name in gene.list){
    # finds pattern 
    if(lazy == FALSE){
      pattern = grep(paste("^", name, sep = ""), 
                     data.table$pretty, 
                     ignore.case = TRUE, 
                     value = TRUE)

    }else if (by.rownames == TRUE){
      pattern = grep(paste("^", name, sep = ""), 
                     rownames(data.table), 
                     ignore.case = TRUE, 
                     value = TRUE)
    }else {
      if (is.null(col) == TRUE){
        stop("No column index given and by.rownames is FALSE")
      }
      pattern = grep(paste("^", name, sep = ""),
                     data.table[ ,col],
                     ignore.case = TRUE,
                     Value = TRUE)
    }
    # if pattern == 0 not found: reutrn message 
    if (length(pattern) == 0 ){
      print(paste(name,": failed to find match", sep = ""))
      next
    }else {
    print(paste("^",name, sep= ""))# prints pattern used
    match <- data.table[pattern,] #matches pattern to dataframe and index
    filtered.data <- rbind(filtered.data, match) # merges to preivous 
    }
  }
  return(filtered.data)
}

