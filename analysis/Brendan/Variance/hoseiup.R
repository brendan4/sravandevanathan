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
    if(log(sd(values)+0.0001) >= -7){
      unexpressed[counter] =  x
      counter = counter+1
    }
  }
  newdata = dataset[unexpressed,]
  return(newdata)
}
#drop_... datasets have unmatched, drop l2_ACAGTG, and l2_ACTTTGA removed
expressed.genes = remove.unexpressed(drop_GeneAbundance)
expressed.trans = remove.unexpressed(drop_t_data)

