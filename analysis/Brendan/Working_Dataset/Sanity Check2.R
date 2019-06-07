library(data.table)
setwd("C:/Users/brendan/Documents/sravandevanathan/analysis/Brendan/working_dataset/Incomplete_Ballgown_data")
old_data <- fread("expressed.genes_WO22.tab")
old_data <- as.data.frame(old_data)
rownames(old_data) <- old_data[,1]
old_data <- old_data[,-c(1)]
colnames(old_data) <- substring(colnames(old_data), first = 6)

new.data <- expressed.genes[,which(colnames(expressed.genes) %in% c(colnames(old_data)))]
colnames(old_data) <- paste("old", colnames(old_data))
new.data$Names <- rownames(new.data)
old_data$Names <- rownames(old_data)
merged.table <- merge(old_data, new.data, by = "Names")
merged.table <- merged.table[,-which( colnames(merged.table) %in% c("Names"))]
merged.cor.table <- cor(na.omit(merged.table), method = "spearman")
merged.cor.table <- as.data.frame(merged.cor.table)
old.col <- merged.cor.table[,grep("^old", colnames(merged.cor.table))]
new.row <- old.col[-grep("^old", rownames(old.col)),]
pairs(new.row, labels = colnames(new.row), upper.panel = NULL)

old_data_sub <- old_data[,c(1,2,3,4,5)]
new_data_sub <- expressed.genes[,which(colnames(expressed.genes) %in% c(colnames(old_data_sub)))]
colnames(old_data_sub) <- paste("old", colnames(old_data_sub))
pairs(as.data.frame(merged.cor.table), labels = colnames(as.data.frame(merged.cor.table)))





new_data_sub$Names <- rownames(new_data_sub)
old_data_sub$Names <- rownames(old_data_sub)
merged.table <- merge(old_data_sub, new_data_sub, by = "Names")
merged.table <- merged.table[,-which( colnames(merged.table) %in% c("Names"))]
                             
pairs(merged.table, labels = colnames(merged.table),upper.panel = NULL)

library(ggplots2)
library(GGally)

merged.cor.table <- cor(na.omit(merged.table), method ="spearman")

pairs(merged.cor.table,labels = colnames(merged.cor.table),upper.panel= NULL)
ggpairs(as.data.frame(merged.table))

