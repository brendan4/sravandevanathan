library(data.table)
setwd("C:/Users/brendan/Documents/sravandevanathan/analysis/Brendan/01_Merged/Old Ballgown")
old_data <- fread("gene_abundance_merged.tab") 
old_data <- as.data.frame(old_data) # changes to data.frame 
old_data <- old_data[,-which ( colnames(old_data) %in% c("V1"))] # drops V1 col 
rownames(old_data) <- old_data[,1] # sets names col as rownames 
old_data <- old_data[,-c(1)] # removes name col
colnames(old_data) <- substring(colnames(old_data), first = 6) # removes FPKM 
old_data <- old_data[,-which ( colnames(old_data) %in% c("L2_ACAGTG", "L6_unmatched","L2_ACTTTGA"))] #removes sample new (and old): optional

#subsets expressed genes
new.data <- expressed.genes[,which(colnames(expressed.genes) %in% c(colnames(old_data)))]
colnames(old_data) <- paste("old", colnames(old_data))
# rownames as col
new.data$Names <- rownames(new.data) 
old_data$Names <- rownames(old_data)
#merging 
merged.table <- merge(old_data, new.data, by = "Names")
merged.table <- merged.table[,-which( colnames(merged.table) %in% c("Names"))]
#cor.table 
merged.cor.table <- cor(na.omit(merged.table), method = "spearman")
merged.cor.table <- as.data.frame(merged.cor.table)
#filtering out unneeded cor values 
old.col <- merged.cor.table[,grep("^old", colnames(merged.cor.table))]
new.row <- old.col[-grep("^old", rownames(old.col)),]

heatmap(as.matrix(new.row))
#heatmap of sample similarities
library(pheatmap)
pheatmap(new.row, cluster_rows = F, cluster_cols = F) 

library(corrplot)
library(grDevices)
corrplot(as.matrix(new.row), cl.lim=c(.85,1), 
         method = "shade", 
         is.cor =F,
         type = "lower",  
         order = "hclust", 
         addrect = 2)

#simple pairs plot ggplot 
library(ggplot2)
library(GGally)
ggpairs(merged.table[,c(3,4,24,25)])
pairs(merged.table[,c(2,3,4,23,24,25)], upper.panel = NULL)

# remove from workspace
rm(old_data,old.col,new.row,merged.cor.table,merged.table,new.data)
