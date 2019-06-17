
library(DBA)

### Gene abundance merge 
main.table <- mergeTables(wd = "C:/Users/brendan/Documents/sravandevanathan/ballgown",
                          commonName = "gene_abundance.tab", 
                          colsToMerge = c(2,5,6,8))

main.table <- merge.cleanup(main.table, cor.table = TRUE)

#removal of L2_ACAGTG
main.table <- main.table[,-which(colnames(main.table) %in% c("FPKM L2_ACAGTG"))]

write.table(main.table, "gene_abundance_merged.tab")
GeneAbundance <- read.table("gene_abundance_merged.tab")

### Transripts merge
main.table <- mergeTables(wd = "C:/Users/brendan/Documents/sravandevanathan/ballgown",
                          commonName = "t_data.ctab", 
                          colsToMerge = c(4,5,6,10,12))

main.table <- merge.cleanup(main.table, cor.table = TRUE)

#removal of L2_ACAGTG
main.table <- main.table[,-which(colnames(main.table) %in% c("FPKM L2_ACAGTG"))]

write.table(main.table,gzfile("transrcipts_merged.ctab.gz")) # writes the table as a .ctab.gz
Transcripts <- read.table(gzfile("transrcipts_merged.ctab.gz")) # will read the table 

rm(main.table)