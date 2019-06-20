#' Pretty gene names
#' 
#' Takes a data.frame with genes append to some unique identifer. EX: RPL11_21321_342342
#' returns RPL11
#' 
#' @param data.set a data.frame with genes to be made pretty
#' @param as.rom.names logical value indicating if the pretty col should turned into row names. 
#' Warning will be triggered if duplicates found and will not execute rownames switch. 
#' @param remove.dups logical statment indicating duplicated rows in pretty col should be removed. 
#' Will always tigger a warning as information across the row is lost. Primarly used to set row names if dups found.
#' @param col an argument only used if the rownames are not where gene names are stored. Should be an integer representing
#' which col should be used. 
#' 
#' @return returns a data.frame with a pretty col containing gene names or a data.frame with pretty rownames. 
#' Rownames will only be added if no duplicates found or remove duplicates specified.
#' @import stringr
#' 
pretty.gene.name <- function(data.set, as.row.names = FALSE, remove.dups = FALSE, col = NULL){
  data.set$pretty <- 0
  if (is.null(col) == TRUE){
    for(row in 1:nrow(data.set)){
      new.name <- str_extract(rownames(data.set[row,]), "^[-a-zA-Z0-9]*")
      data.set$pretty[row] <- new.name 
    }
  }else{
    for(row in 1:nrow(data.set)){
      new.name <- str_extract(data.set[row,col], "^[-a-zA-Z0-9]*")
      data.set$pretty[row] <- new.name 
  }
  }
  if (remove.dups == TRUE){
    warning("Removing duplicates: information will be lost!")
    data.set <- data.set[-which(duplicated(data.set$pretty)),]
  }
  if(as.row.names == TRUE){
    if (sum(duplicated(data.set$pretty)) > 0){
      warning("Duplicated values found")
    } else{
      rownames(data.set) <- data.set$pretty
      data.set <- data.set[,-which(colnames(data.set) %in% c("pretty"))]
    }
  }
  return(data.set)
}



