
library(DBA)

#gene abundance merge
main.table <- mergeTables(wd = "C:/Users/brendan/Documents/sravandevanathan/ballgown",
                          commonName = "gene_abundance.tab", 
                          colsToMerge = c(2,5,6,8))

GeneAbundance <- merge.cleanup(main.table, cor.table = TRUE)

