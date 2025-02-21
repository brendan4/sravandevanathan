library(data.table)
library(tidyr)
library(devtools)
library(DBA)
library(pheatmap)
library(ape)

### Gene abundance merge 
main.table <- mergeTables(wd = "~/sravandevanathan/ballgown_new annotation",
                          commonName = "gene_abundance.tab", 
                          colsToMerge = c(2,5,6,8))

setwd("~/sravandevanathan/analysis/Brendan") 

data <- merge.cleanup(main.table,
                      boxplot = F,
                      cor.table = TRUE,
                      tidy.colnames = TRUE,
                      remove.NA = FALSE)

cor.table <- data[[1]]
main.table <- data[[2]]

#removal of unmatched 
main.table <- main.table[,-which(colnames(main.table) %in% 
                                   c("L3_unmatched", "L6_unmatched"))]

#correlation plots
colnames(main.table)[1:8] <- substring(colnames(main.table)[1:8], first = 8)
cor.plots(main.table, heatmap = TRUE, phylo = FALSE)

main.table <- main.table[,-which(colnames(main.table) %in% 
                                   c("L2_ACAGTG", "LIB4-98412318", "LIB7-98397315",
                                     "LIB8-98401312", "LIB2-98397314", 
                                     "L3_ATCACG", "L3_TTAGGC", "L6_TTAGGC"))]

cor.plots(main.table, heatmap = TRUE, phylo = FALSE)

write.table(main.table, "gene_abundance_merged.tab")
GeneAbundance <- read.table("gene_abundance_merged.tab")

### Transripts merge
main.table <- mergeTables(wd = "~/sravandevanathan/ballgown_new annotation",
                          commonName = "t_data.ctab", 
                          colsToMerge = c(4,5,6,10,12))

setwd("~/sravandevanathan/analysis/Brendan")

data <- merge.cleanup(main.table,
                      cor.table = TRUE,
                      tidy.colnames = TRUE)

cor.table <- data[[1]]
main.table <- data[[2]]

#removal of unmatched 
main.table <- main.table[,-which(colnames(main.table) %in% 
                                   c("L3_unmatched", "L6_unmatched"))]

#correlation plots
colnames(main.table)[1:8] <- substring(colnames(main.table)[1:8], first = 8) # cleans July lib names
cor.plots(main.table, heatmap = TRUE, phylo = FALSE)

#removal of L2_ACAGTG
main.table <- main.table[,-which(colnames(main.table) %in% 
                                   c("L2_ACAGTG", "LIB4-98412318", "LIB7-98397315",
                                     "LIB8-98401312", "LIB2-98397314", 
                                     "L3_ATCACG", "L3_TTAGGC", "L6_TTAGGC"))]

cor.plots(main.table, heatmap = TRUE, phylo = FALSE)

write.table(main.table,gzfile("transrcipts_merged.ctab.gz")) # writes the table as a .ctab.gz
Transcripts <- read.table(gzfile("transrcipts_merged.ctab.gz")) # will read the table 

#variation 
plot.var(GeneAbundance)
plot.var(Transcripts, title.name = "Transcripts")
abline(v = -7, col = "red")

#remove unexpressed
expressed.genes = remove.unexpressed(GeneAbundance, cutoff = -7) # cutoff at -7
expressed.trans = remove.unexpressed(Transcripts, cutoff = -7) # cutoff at -7 
expressed.genes <- na.omit(expressed.genes)
expressed.trans <- na.omit(expressed.trans)

# saving new data
setwd("C:/Users/brendan/Documents/sravandevanathan/")
write.table(expressed.genes, "expressed.genes.GEN.tab")
write.table(expressed.trans, "expressed.trans.GEN.tab")

clean.environment()

# hierarchical clustering of samples and heatmap of sample similarities: zoom for heatmap
cor.plots(expressed.genes, heatmap = TRUE, phylo = TRUE)
cor.plots(expressed.trans, heatmap = TRUE, phylo = TRUE)

#PCA and MDS
PCA(expressed.genes.GEN, scaled = FALSE, PCA.Genes = FALSE)
PCA(expressed.trans)
MDS(expressed.genes, scaled = FALSE)
MDS(expressed.trans)

data("expressed.genes")