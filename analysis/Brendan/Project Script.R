library(data.table)
library(tidyr)
library(devtools)
library(DBA)
library(pheatmap)
library(ape)

### Gene abundance merge 
main.table <- mergeTables(wd = "~/sravandevanathan/ballgown",
                          commonName = "gene_abundance.tab", 
                          colsToMerge = c(2,5,6,8))

data <- merge.cleanup(main.table,
                      boxplot = TRUE,
                      cor.table = TRUE,
                      tidy.colnames = TRUE,
                      remove.NA = FALSE)

cor.table <- data[[1]]
main.table <- data[[2]]

#removal of L2_ACAGTG
main.table <- main.table[,-which(colnames(main.table) %in% 
                                   c("FPKM L2_ACAGTG","FPKM.L3_unmatched","FPKM.L6_unmatched"))]
write.table(main.table, "gene_abundance_merged.tab")
GeneAbundance <- read.table("gene_abundance_merged.tab")

### Transripts merge
main.table <- mergeTables(wd = "~/sravandevanathan/ballgown",
                          commonName = "t_data.ctab", 
                          colsToMerge = c(4,5,6,10,12))

data <- merge.cleanup(main.table,
                      cor.table = TRUE,
                      tidy.colnames = TRUE)

cor.table <- data[[1]]
main.table <- data[[2]]

#removal of L2_ACAGTG
main.table <- main.table[,-which(colnames(main.table) %in% 
                                   c("FPKM L2_ACAGTG","FPKM.L3_unmatched","FPKM.L6_unmatched"))]
write.table(main.table,gzfile("transrcipts_merged.ctab.gz")) # writes the table as a .ctab.gz
Transcripts <- read.table(gzfile("transrcipts_merged.ctab.gz")) # will read the table 

#variation 
plot.var(GeneAbundance)
plot.var(Transcripts)
abline(a = 0, b = 0, v = -7, col = "red")

#remove unexpressed
expressed.genes = remove.unexpressed(GeneAbundance, cutoff = -7) # cutoff at -7
expressed.trans = remove.unexpressed(Transcripts, cutoff = -8) # cutoff at -8 

# saving new data
setwd("C:/Users/brendan/Documents/sravandevanathan/")
write.table(expressed.genes, "expressed.genes.tab")
write.table(expressed.trans, "expressed.trans.tab")

clean.environment()

# hierarchical clustering of samples and heatmap of sample similarities: zoom for heatmap
cor.plots(expressed.genes, heatmap = TRUE, phylo = TRUE)
cor.plots(expressed.trans, heatmap = TRUE, phylo = TRUE)

#PCA and MDS
PCA(expressed.genes, scaled = FALSE, PCA.Genes = FALSE)
PCA(expressed.trans)MDS(expressed.genes, scaled = FALSE)

MDS(expressed.trans)
