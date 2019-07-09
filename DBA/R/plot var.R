plot.var = function(dataset, zero.offset = 0.0001){
  var <- apply(dataset,1, sd)
  plot(density(log(var + zero.offset)), xlim=c(-15,10))
}
