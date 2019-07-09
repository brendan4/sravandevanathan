library(RColorBrewer)

#annotation setup
data("full.pheno.table")
pheno.table <- full.pheno.table[,c(2,1)]
rownames(pheno.table) <- pheno.table[,2]
pheno.table[,2] <- full.pheno.table[,3]
colnames(pheno.table) <- c("Pheno", "Replicates")

# Specify colors for annnotation: Phenotype
Pheno = c("red", "yellow", "grey", "black")
names(Pheno) = c("SS", "S", "C", "W")

# Specify colors for annnotation: Replicates
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
Reps = c(col=sample(col_vector, 15))
names(Reps) = full.pheno.table[!(duplicated(full.pheno.table$Replicates)), "Replicates"]

# preping data
ann_colors = list(Pheno = Pheno, Replicates = Reps)

#####heatmap corr between samples

cor.plots(expressed.genes, 
          method = "spearman", 
          annotation = pheno.table, 
          colors = ann_colors)

#### Individauls PCA
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
PCA(simple.names, pheno = pheno, label.size = 4, pca.dim = c(2,3), scaled = FALSE, color.option = 1)

#mystery samples: found in all three lanes but not in pheno excel 
"L3_TTAGGC" %in% colnames(expressed.genes)
grep("TTAGGC", colnames(expressed.genes), value = T)

#### PCA with pheno data only 
PCA(expressed.genes, pheno = pheno, label.size = 3)

colors <- c("green", "red", "orange", "blue")
names(colors) <- c(pheno$pheno[!duplicated(pheno$pheno)])