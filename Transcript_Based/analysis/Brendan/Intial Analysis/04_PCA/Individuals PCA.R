data("full.pheno")
data("pheno")

#changes samples names to individauls names 
change.names <- function(data.set, full.pheno, pheno){
  pheno<- t(pheno)
  for(person in 1:length(full.pheno$Wildtype)){
    name <- names(full.pheno$Wildtype)[person]
    samples <- full.pheno$Wildtype[[person]]
    curr.names <- colnames(data.set)
    curr.names[which(curr.names %in% samples)] <- name
    pheno[1,which(pheno[1,] %in% samples)] <- name
    colnames(data.set) <- curr.names
  }
  for (subcat in 1:length(full.pheno$Mutant)){
    for (person in 1:length(full.pheno$Mutant[[subcat]])){
      name <- names(full.pheno$Mutant[[subcat]])[person]
      samples <- full.pheno$Mutant[[subcat]][[person]]
      curr.names <- colnames(data.set)
      curr.names[which(curr.names %in% samples)] <- name
      pheno[1,which(pheno[1,] %in% samples)] <- name
      colnames(data.set) <- curr.names
    }
  }
  pheno <- t(pheno)
  results <- list(data.set, pheno)
  return(results)
}

simple <- change.names(expressed.genes, full.pheno, pheno)
simple.pheno <- simple[[2]]
simple.names <- simple[[1]]

#PCA 
PCA(simple.names, pheno = pheno, label.size = 4, pca.dim = c(1,2), scaled = FALSE)
#mystery samples
"L3_TTAGGC" %in% colnames(expressed.genes)
#found in all three lanes but not in the excel 
grep("TTAGGC", colnames(expressed.genes), value = T)



