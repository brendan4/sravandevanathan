#filers data set with genelist using rownames 
filter.genes <- function(data.table, gene.list){
  filtered.data <- data.frame()
  for (name in gene.list){
    # finds pattern 
    pattern = grep(paste("^", name, sep = ""), 
                   rownames(data.table), 
                   ignore.case = TRUE, 
                   value = TRUE)
    # if pattern == 0 not found: reutrn message 
    if (pattern == 0 ){
      print(paste(pattern,": failed to find match", sep = ""))
      next
    } else {
    print(paste("^",name, sep= ""))# prints pattern used
    match <- data.table[pattern,] #matches pattern to dataframe and index
    filtered.data <- rbind(filtered.data, match) # merges to preivous 
    }
  }
  return(filtered.data)
}

RPL <- filter.genes(expressed.genes, "RPL")
