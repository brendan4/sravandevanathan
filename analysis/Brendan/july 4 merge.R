

july <- mergeTables(wd = "~/sravandevanathan/ballgown_july_4", 
            commonName = "gene_abundance.tab.gz",
            colsToMerge = c(2,5,6,8))

summary(july)
merge.cleanup(july, boxplot = TRUE, cor.table = TRUE)
