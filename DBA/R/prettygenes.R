pretty.gene.name <- function(data.set){
  
  data.set$pretty <- 0
  for(row in 1:nrow(data.set)){
    new.name <- str_extract(rownames(data.set[row,]), "^[-a-zA-Z0-9]*")
    data.set$pretty[row] <- new.name 
    
  }
  
  return(data.set)
}



