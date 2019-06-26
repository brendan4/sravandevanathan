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
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
Reps = c(col=sample(col_vector, 16))
names(Reps) = full.pheno.table[!(duplicated(full.pheno.table$Replicates)), "Replicates"]

# preping data
ann_colors = list(Pheno = Pheno, Replicates = Reps)

#heatmap corr between samples
cor.plots(expressed.genes, method = "spearman", annotation = pheno.table, colors = ann_colors)
