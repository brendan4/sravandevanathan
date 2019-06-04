remove.unexpressed = function(dataset){
  unexpressed = c()
  counter = 1
  for (x in 1:nrow(dataset)){
    values = c()
    for (y in 1:ncol(dataset)-1){
      if(is.na(dataset[x,y+1]) == T){
        values[y] = 0.0
      }else{
        values[y] = dataset[x,y+1]
      }
    }
    if(log(sd(values)+0.0001) >= -8){
      unexpressed[counter] =  x
      counter = counter+1
    }
  }
  newdata = dataset[unexpressed,]
  return(newdata)
}
#drop_... datasets have unmatched, drop l2_ACAGTG removed
expressed.genesN = remove.unexpressed(GeneAbundance)
expressed.transN = remove.unexpressed(Transcripts)


# Anaylsis after first filter: drop L6_unmatched + L3_ATCACG
boxplot(log2(na.omit(expressed.transN)+ 0.1), 
        names=colnames(expressed.transN), las=2, ylab="log2(FPKM)", 
        main="Distribution of FPKMs for all libraries")
summary(expressed.transN)

