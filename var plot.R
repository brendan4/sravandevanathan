#drop l2_ACAGTG and l2_ACTTTGA
drop_GeneAbundance <- GeneAbundance[c(-2,-4)]
drop_t_data <- t_data[c(-2,-4)]

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
    variances[x] =  log(sd(values))
  }
  plot(density(variances),xlim=c(-10,20))
}

#call to var func
plot.var(drop_GeneAbundance)
plot.var(drop_t_data)
