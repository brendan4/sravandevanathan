
july <- mergeTables(wd = "~/sravandevanathan/ballgown_july_4", 
            commonName = "gene_abundance.tab.gz",
            colsToMerge = c(2,5,6,8))

july <- as.data.frame(july)
summary(july)

items <- merge.cleanup(july, boxplot = TRUE, cor.table = TRUE, tidy.colnames = TRUE)
cor.table <- items[[1]]
july <- items[[2]]

# 121317_LIB4 most disimilar = II.8 in phenotable
data("full.pheno.table")
II.8 <- full.pheno.table[grep("II.eight", full.pheno.table$Replicates),1]
II.8 <- expressed.genes[grep(II.8,colnames(expressed.genes))]
II.8.data <- merge(II.8, july[ ,"121317_LIB4-98412318",drop = FALSE], by = "row.names", all.x=TRUE)
cor(na.omit(II.8.data[,-c(1)]))

mean(na.omit(july$`121317_LIB4-98412318`))
sd(na.omit(july$`121317_LIB4-98412318`))

#removal of "121317_LIB4-98412318"
july <- july[,-which(colnames(july) %in% "121317_LIB4-98412318")]

#cleaning 
colnames(july) <- substring(colnames(july), 
                                first = 8)

#variation 
plot.var(july)
abline(a = 0, b = 0, v = -6, col = "red")

#variation filtering 
july <- remove.unexpressed(july, cutoff = -6)

#correlation plots 
cor.plots(july, heatmap = TRUE, phylo = TRUE, log = TRUE)

#PCA plots
PCA(july, scaled = FALSE, PCA.Genes = FALSE, label.size = 4)
MDS(july, scaled = TRUE)

cor(na.omit(july), method = "spearman")

#merging with expressed.genes data
new <- merge(expressed.genes, july, by = "row.names", all.x = TRUE, all.y = TRUE)
row.names(new) <- new[,1]
new <- new[,-c(1)]
cor.plots(na.omit(new), heatmap = TRUE)

#throwing out the samples LIB8, LIB2, LIB7
new <- new[,-which(colnames(new) %in% c("LIB8-98401312", 'LIB2-98397314', 'LIB7-98397315'))]
july <- july[,-which(colnames(july) %in% c("LIB8-98401312", 'LIB2-98397314', 'LIB7-98397315'))]
new <- new[, -which(colnames(new) %in% c('means', 'pretty'))]

#cor plot after filtering 
cor.plots(na.omit(new), heatmap = TRUE)
new.save <- new 

#not finished
#changing indivdual names 
for (i in 1:length(new)){
  individaul <- all.pheno.data[which(all.pheno.data$`colnames(expressed.genes)` %in% colnames(new)[i]),3]
  colnames(new)[i] <- individaul
}

PCA(new, pca.dim = c(1,2), color.option = 1, pheno = all.pheno.data, label.size = 4)

all.pheno.data[,c(1,2)]
