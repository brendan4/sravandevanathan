
merge.cleanup <- function(data.set, 
                          boxplot = FALSE, 
                          cor.table = FALSE, 
                          remove.NA = FALSE, 
                          tidy.colnames = FALSE){
  
  data.set <- as.data.frame(data.set) # forms data.table to data.frame
  rownames(data.set) <- data.set[,1] #gene namess to row names
  data.set <- data.set[,-c(1)] # removes gene names col
  
  if (boxplot == TRUE){
  boxplot(log2(na.omit(data.set)+.1), 
          names=colnames(data.set), las=2, ylab="log2(FPKM)", 
          main="Distribution of FPKMs for all libraries")
  }
  
  if (cor.table == TRUE){
    temp <- na.omit(data.set)
    cor.gene <- cor(temp, method = "spearman")
    print(cor.gene)
  }
  
  if (remove.NA == TRUE){
    data.set <- na.omit(data.set)
  }
  
  if(tidy.colnames == TRUE){
    colnames(data.set) <- substring(colnames(data.set), 
                                    first = 6) # removes FPKM. from row names 
  }
  return(data.set)
}

