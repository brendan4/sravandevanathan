#opening GTEX median expression for cell types: data loadings
library(data.table)
library(ggplot2)

setwd("~/sravandevanathan/analysis/Brendan/Cellular Heterogenity/xCell/Sanity/Datasets")
GTEX <- fread("GTEx_gene_median_tpm.gct.gz")
GTEX <- as.data.frame(GTEX)
GTEX <- GTEX[,which(colnames(GTEX) %in% c("Description","Whole Blood"))]

means <- apply(expressed.genes, 1, mean)
means <- na.omit(means)

high.means <- means[which(means > 800)]
high <- filter.genes(expressed.genes, names(high.means))
high.GTEX <- GTEX[GTEX$`Whole Blood` > 2000,]

high <- pretty.gene.name(high, as.row.names = TRUE)

high.GTEX$Description[which(!high.GTEX$Description %in% high$pretty)]
high$pretty[which(!high$pretty %in% high.GTEX$Description)]
