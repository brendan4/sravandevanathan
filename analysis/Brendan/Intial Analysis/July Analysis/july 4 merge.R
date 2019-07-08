

july <- mergeTables(wd = "~/sravandevanathan/ballgown_july_4", 
            commonName = "gene_abundance.tab.gz",
            colsToMerge = c(2,5,6,8))

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

#correlation plots 
cor.plots(july, heatmap = TRUE, phylo = TRUE, log = TRUE)

#PCA plots
PCA(july, scaled = FALSE, PCA.Genes = FALSE, label.size = 4)
MDS(july, scaled = TRUE)


#dropping L6_TTAGGC and L3_TTAGGC
df <- df[-which(df$id %in% c("L6_TTAGGC")),]
df <- df[-which(df$person  == 0), ]