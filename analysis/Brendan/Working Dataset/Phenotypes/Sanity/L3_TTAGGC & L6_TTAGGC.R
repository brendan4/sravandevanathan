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

data("full.pheno.table")
expressed.genes <- expressed.genes[, - which(col.names(expressed.genes) %in% c("L3_TTAGGC", "L6_TTAGGC")]