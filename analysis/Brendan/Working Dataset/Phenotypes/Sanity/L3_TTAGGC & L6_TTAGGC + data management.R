# possible mistake in phenotype table indivduals correlate however they do not appear to have same sex 
sub <- expressed.genes[,which(colnames(expressed.genes) %in% c("L6_TTAGGC", "L3_ACTTGA"))]
sex.genes <- c("XIST", "USP9Y", "UTY", "RPS4Y1", "TSIX", "PRKY", 'DDX3Y', 'RPS4Y1')
sub.filt <- filter.genes(sub, sex.genes)
sub.hem <- filter.genes(sub,"hb")
setwd("C:/Users/brendan/Documents/sravandevanathan/analysis/Brendan/Working Dataset/02_Filtered Data")
filt.names <- read.delim("sex_related.tab", sep= "")
sub.filt <- filter.genes(sub, filt.names$geneNames)

cor(na.omit(sub), method = "spearman")

#sex profile of L3_TTAGGC
sex.genes <- c("XIST", "USP9Y", "UTY", "RPS4Y1", "TSIX", "PRKY", 'DDX3Y', 'RPS4Y1')
TTAGGC <- expressed.genes[,grep("TTAGGC", colnames(expressed.genes))]
sub.filt <- filter.genes(TTAGGC, sex.genes)

#corelation to other TTAGGC
cor(na.omit(TTAGGC), method = "spearman")

# corelation to L3_GATCAG individual II.6 (had same barcode in earlier hiseq run)
sub <- expressed.genes[,which(colnames(expressed.genes) %in% c("L3_TTAGGC", "L3_GATCAG"))]
sub.filt <- filter.genes(sub, sex.genes)
cor(na.omit(sub), method = "spearman")


#droppping L3_TTAGGC and L6_TTAGGC from all data
expressed.genes <- expressed.genes[, - which(colnames(expressed.genes) %in% c("L3_TTAGGC", "L6_TTAGGC"))]
expressed.trans <- expressed.trans[, -which(colnames(expressed.trans) %in% c("L3_TTAGGC", "L6_TTAGGC"))]

#redoing all cutoff
plot.var(expressed.genes)
abline(a = 0, b = 0, v = -7, col = "red")
expressed.genes <- remove.unexpressed(expressed.genes, cutoff = -7)
plot.var(expressed.trans)
abline(a = 0, b = 0, v = -8, col = "red")
expressed.trans <- remove.unexpressed(expressed.trans, cutoff = -8)

data("pheno.basic")
data("pheno.colors")
data("full.pheno")
data("pheno")
data("full.pheno.table")

pheno.basic <- pheno.basic[-which(pheno.basic$`colnames(expressed.genes)` %in% c("L3_TTAGGC", "L6_TTAGGC")), ]
full.pheno.table <- full.pheno.table[-which(full.pheno.table$`colnames(expressed.genes)` %in% c("L3_TTAGGC", "L6_TTAGGC")), ]
pheno <- pheno[-which(pheno$`colnames(expressed.genes)` %in% c("L3_TTAGGC", "L6_TTAGGC")), ]
full.pheno$Mutant$Carrier$I.one[which(full.pheno$Mutant$Carrier$I.one %in% c("L3_TTAGGC", "L6_TTAGGC"))] <- 0 
pheno.colors <- pheno.colors[-which(names(pheno.colors) %in% c("L3_TTAGGC", "L6_TTAGGC"))]

usethis::use_data(expressed.trans, DBA, overwrite = T)
usethis::use_data(expressed.genes, DBA, overwrite = T)
usethis::use_data(pheno.colors, DBA, overwrite = TRUE)
usethis::use_data(pheno, DBA, overwrite = TRUE)
usethis::use_data(pheno.basic, DBA,overwrite = TRUE)
usethis::use_data(full.pheno, DBA, overwrite = TRUE)
usethis::use_data(full.pheno.table, DBA, overwrite = TRUE)
