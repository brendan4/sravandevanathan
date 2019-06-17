plot.var = function(dataset, zero.offset = 0.0001){
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