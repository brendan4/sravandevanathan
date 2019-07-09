remove.unexpressed = function(dataset, cutoff, zero.offset = 0.0001){
  dataset[is.na(dataset)] <- 0
  var <- apply(dataset,1, sd)
  var <- log(var + zero.offset)
  idx <- which(var >= cutoff)
  newdata = dataset[idx,]
  return(newdata)
}
