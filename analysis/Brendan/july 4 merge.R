

july <- mergeTables(wd = "~/sravandevanathan/ballgown_july_4", 
            commonName = "gene_abundance.tab.gz",
            colsToMerge = c(2,5,6,8))

summary(july)
items <- merge.cleanup(july, boxplot = TRUE, cor.table = TRUE, tidy.colnames = TRUE)
cor.table <- items[[1]]
test <- items[[2]]
