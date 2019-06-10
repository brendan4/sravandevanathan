summary(GeneAbundance)
summary(Transcripts)

cor.gene <- cor(GeneAbundance, method = "spearman")
print("Gene: spearman correlation matrix")
print(cor.gene)
cor.trans <- cor(Transcripts, method = "spearman")
print("Transcripts: spearman correlation matrix")
print(cor.gene)

#filter out both unmatched and L3_ATCACG
print("Flitering out FPKM.L3_unmatched, FPKM.L3_ATCACG, FPKM.L6_unmatched")
GeneAbundance <- GeneAbundance[,-which(colnames(GeneAbundance) %in% c("FPKM.L3_unmatched", "FPKM.L3_ATCACG", "FPKM.L6_unmatched"))]
Transcripts <- Transcripts[,-which(colnames(Transcripts) %in% c("FPKM.L3_unmatched", "FPKM.L3_ATCACG","FPKM.L6_unmatched"))]

plot.var = function(dataset){
  variances = c()
  for (x in 1:nrow(dataset)){
    values = c()
    for (y in 1:ncol(dataset)-1){
      if(is.na(dataset[x,y+1]) == T){
        values[y] = 0.0
      }else{
        values[y] = dataset[x,y+1]
      }
    }
    variances[x] =  log(sd(values)+0.0001)
  }
  plot(density(variances),xlim=c(-15,10))
}

plot.var(na.omit(GeneAbundance))
plot.var(na.omit(Transcripts))
abline(a = 0, b = 0, v = -7, col = "red")

remove.unexpressed = function(dataset, cutoff){
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
    if(log(sd(values)+0.0001) >= cutoff){
      unexpressed[counter] =  x
      counter = counter+1
    }
  }
  newdata = dataset[unexpressed,]
  return(newdata)
}

print("filtering datasets")
expressed.genesN = remove.unexpressed(GeneAbundance, cutoff = -7) # cutoff at -7
expressed.transN = remove.unexpressed(Transcripts, cutoff = -8) # cutoff at -8 
print("remove unfiltered data")
rm(GeneAbundance)
rm(Transcripts)

# saving new data
setwd("C:/Users/brendan/Documents/sravandevanathan/")
write.table(expressed.genes, "expressed.genes.tab")
write.table(expressed.trans, "expressed.trans.tab")