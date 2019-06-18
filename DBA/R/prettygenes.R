pretty.gene.name <- function(data.set, as.row.names = FALSE, remove.dups = FALSE){
  
  data.set$pretty <- 0
  for(row in 1:nrow(data.set)){
    new.name <- str_extract(rownames(data.set[row,]), "^[-a-zA-Z0-9]*")
    data.set$pretty[row] <- new.name 
    
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



