
merge.cleanup <- function(data.set, boxplot = FALSE, cor.table = FALSE){
  data.set <- as.data.frame(data.set) # forms data.table to data.frame
  rownames(data.set) <- data.set[,1] #gene namess to row names
  data.set <- data.set[,-c(1)] # removes gene names col
  
  if (boxplot == TRUE){
  boxplot(log2(na.omit(main.table)+.1), 
          names=colnames(main.table), las=2, ylab="log2(FPKM)", 
          main="Distribution of FPKMs for all libraries")
  }
  
  if (cor.table == TRUE){
    cor.gene <- cor(data.set, method = "spearman")
    print(cor.table)
  }
  return(data.set)
}

