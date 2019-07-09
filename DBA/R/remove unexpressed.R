remove.unexpressed = function(dataset, cutoff, zero.offset = 0.0001){
  dataset.na <- dataset
  dataset.na[is.na(dataset)] <- 0
  var <- apply(dataset.na,1, sd)
  var <- log(var + zero.offset)
  idx <- which(var >= cutoff)
  newdata = dataset[idx,]
  return(newdata)
}
