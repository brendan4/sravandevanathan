#drop l2_ACAGTG and l2_ACTTTGA and unmatched 
drop_GeneAbundance <- GeneAbundance[c(-2,-4,-25)]
drop_t_data <- t_data[c(-2,-4,-25)]
rownames(drop_GeneAbundance) <- drop_GeneAbundance[,1]
rownames(drop_t_data) <- drop_t_data[,1]
drop_GeneAbundance <- drop_GeneAbundance[-1]
drop_t_data <- drop_t_data[-1]

#var function plot 
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

#call to var func
plot.var(GeneAbundance)
plot.var(Transcripts)
abline(a = 0, b = 0, v = -7, col = "red")
