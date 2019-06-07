#computes correlation and returns a matrix
cor.table <- function(data.table){
  cor.matrix <- data.frame()
  
  # generate an empty grid 
  for(i in 1:length(data.table)){
    for(j in 1:length(data.table)){
      cor.matrix[i,j] <- 0
    }
  }
  
  
  rownames(cor.matrix) <- sort(colnames(data.table))
  colnames(cor.matrix) <- sort(colnames(data.table))
  
  #computation of correlation 
  for (i in 1:length(data.table)){
    for (j in 1:length(data.table)){
      # skipp if equal 
      if (i == j){
        cor.matrix[i,j] <- 1
      }else{
        x = data.table[,i]
        y = data.table[,j]
        rs <- cor(x,y, method="spearman")
        cor.matrix[i,j] <- rs 
      }
      
    }
  
  }
  return(cor.matrix)
  
}


cor.matrix  <- cor.table(na.omit(expressed.genes)) # testing with expressed.genes data 

#testing or matrix 
test <- na.omit(expressed.genes)
x = test$L2_ACTTGA
y = test$L2_ATCACG
rs <- cor(x,y)^2

rs == cor.matrix["L2_ACTTGA","L2_ATCACG"]
rs == cor.matrix["L2_ATCACG", "L2_ACTTGA"]
cor.matrix["L2_ATCACG", "L2_ACTTGA"] == cor.matrix["L2_ACTTGA","L2_ATCACG"]

rm(rs, x, y, test) # remove test from enviroment 
rm(cor.matrix) # removes output
rm(cor.table) # removes function 

