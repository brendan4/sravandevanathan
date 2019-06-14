#filers data set with genelist using rownames 
filter.out.genes <- function(data.table, gene.list){
  filtered.data <- data.table
  for (name in gene.list){
    print(name)
    info <- grepl(paste("^",name, sep = ""),
                   rownames(filtered.data))
    
    
    filtered.data <- filtered.data[which(info == F),] #matches pattern to dataframe and index
  
  
  }
  return(filtered.data)
}

