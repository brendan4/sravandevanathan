data("pheno")
data("full.pheno")

sample.names <- names(full.pheno$Wildtype)

for(sample in 1:length(full.pheno$Wildtype)){
  id.holder <- unlist(full.pheno$Wildtype[sample], use.names=FALSE)
  cols <- c(which(colnames(expressed.genes) %in% c(id.holder)))
  print(names(full.pheno$Wildtype[sample]))
  print(colnames(expressed.genes)[cols])
    }

full.pheno$Wildtype$III.one[which(!(full.pheno$Wildtype$III.one %in% colnames(expressed.genes)))]

validation <- function(full.pheno, data.set){
  num.lost <- 0
  lost.values <- c()
  
  #Wildtype
  for (sample in 1:length(full.pheno$Wildtype)){
    id.holder <- unlist(full.pheno$Wildtype[sample], use.names = FALSE)
    lost<- which(!(id.holder %in% colnames(expressed.genes)))
    if(length(lost) == 0){
      print(paste(names(full.pheno$Wildtype)[sample], ":", full.pheno$Wildtype[sample], ": is good"), sep= " ")
    } else{
      num.lost <- num.lost + length(lost)
      lost <- names(full.pheno$Wildtype[sample][c(lost)])
      lost.values <- append(lost.values, lost)
      print(lost)
    }
  }
  
  for (cat in 1:length(full.pheno$Mutant)){
    for (subcat in 1:length(full.pheno[cat])){
      id.holder <- unlist(full.pheno$Mutant[cat][subcat], use.names = FALSE)
      lost <- which
      
    }
  }
  stop()
  print(paste("Total lost samples: ", toString(num.lost)))
  return(lost.values)
}
validation(full.pheno = full.pheno, expressed.genes)
