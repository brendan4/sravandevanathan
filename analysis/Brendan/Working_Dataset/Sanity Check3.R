library(data.table)
setwd("C:/Users/brendan/Documents/sravandevanathan/analysis/Brendan/merged/Old Ballgown")
old_data <- fread("gene_abundance_merged.tab") 
old_data <- as.data.frame(old_data) # changes to data.frame 
old_data <- old_data[,-which ( colnames(old_data) %in% c("V1"))] # drops V1 col 
rownames(old_data) <- old_data[,1] # sets names col as rownames 
old_data <- old_data[,-c(1)] # removes name col
colnames(old_data) <- substring(colnames(old_data), first = 6) # removes FPKM 
old_data <- old_data[,-which ( colnames(old_data) %in% c("L2_ACAGTG", "L6_unmatched"))] #removes sample new (and old): optional

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

# remove from workspace
rm(old_data,old.col,new.row,merged.cor.table,merged.table,new.data)
