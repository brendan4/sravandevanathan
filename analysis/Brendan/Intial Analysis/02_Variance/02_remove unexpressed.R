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
#l2_ACAGTG removed
expressed.genesN = remove.unexpressed(GeneAbundance)
expressed.transN = remove.unexpressed(Transcripts)


# Anaylsis after first filter: drop L6_unmatched + L3_ATCACG
boxplot(log2(na.omit(expressed.genes)+ 0.1), 
        names=colnames(expressed.genes), las=2, ylab="log2(FPKM)", 
        main="Distribution of FPKMs for all libraries")
summary(expressed.genesN)

#both unmatched removed L3_ATCACG removed 
expressed.genesN = remove.unexpressed(GeneAbundance) # cutoff at -7
expressed.transN = remove.unexpressed(Transcripts) # cutoff at -8 
expressed.genes <- expressed.genesN
expressed.trans <- expressed.transN
rm(expressed.genesN)
rm(expressed.transN)
rm(GeneAbundance)
rm(Transcripts)

setwd("C:/Users/brendan/Documents/sravandevanathan/")
write.table(expressed.genes, "expressed.genes.tab")
write.table(expressed.trans, "expressed.trans.tab")

