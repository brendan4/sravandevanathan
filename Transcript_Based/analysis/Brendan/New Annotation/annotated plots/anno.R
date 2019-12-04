library(RColorBrewer)
#annotation setup
data("expressed.genes.GEN")
data("pheno.more.samples")
pheno.table <- pheno.more.samples[,c(2,1)]
rownames(pheno.table) <- pheno.table[,2]
pheno.table[,2] <- pheno.more.samples[,3]
colnames(pheno.table) <- c("Pheno", "Replicates")

# Specify colors for annnotation: Phenotype
Pheno = c("red", "yellow", "grey", "black")
names(Pheno) = c("SS", "S", "C", "W")

# Specify colors for annnotation: Replicates
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
Reps = c(col=sample(col_vector, 16))
names(Reps) = pheno.more.samples[!(duplicated(pheno.more.samples$Replicates)), "Replicates"]

# preping data
ann_colors = list(Pheno = Pheno, Replicates = Reps)

#####heatmap corr between samples

cor.plots(expressed.genes.GEN, 
          method = "spearman", 
          annotation = pheno.table, 
          colors = ann_colors)


#### Individauls PCA

data("pheno.more.samples")
ind <- pheno.more.samples[,"Replicates"]
ind <- ind[!duplicated(ind)]
for(i in 1:length(ind)){
  samples <- pheno.more.samples[which(pheno.more.samples$Replicates %in% ind[i]), "colnames(expressed.genes)"]
  colnames(expressed.genes.GEN)[which(colnames(expressed.genes.GEN) %in% samples)] <- ind[i]
}

#PCA replicates 
look <- PCA(expressed.genes.GEN, 
    pheno = pheno.more.samples[,1:2], 
    label.size = 4, 
    pca.dim = c(1,2), 
    scaled = T, 
    color.option = 1,
    PCA.Genes = TRUE, title = "Gene Level - log")

PCA.look <- filter.genes(expressed.genes.GEN, gene.list = look)
expressed.genes.GEN <- filter.out.genes(expressed.genes.GEN, gene.list = "^RN7")
 