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
individuals <- full.pheno.table[duplicated(full.pheno.table$Replicates), "Replicates"]
cor.data <- list()
mean.data <- c()

for(i in 1:length(individuals)){
  sub <- full.pheno.table[which(full.pheno.table$Replicates == individuals[i]),]
  sub.data <- expressed.genes[,which(colnames(expressed.genes) %in% sub[,1])]
  cor.table <- cor(na.omit(sub.data), method = "spearman")
  cor.data[[i]] <- cor.table
  names(cor.data)[i] <- individuals[i]
  mean.data[i] <- mean(cor.table)
  names(mean.data)[i] <- individuals[i]
}

# least cor between samples
cor.data[[5]]
names(cor.data)[5]



