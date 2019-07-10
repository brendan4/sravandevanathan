#opening GTEX median expression for cell types: data loadings
library(data.table)
library(ggplot2)

setwd("~/sravandevanathan/analysis/Brendan/Cellular Heterogenity/xCell/Sanity/Datasets")
GTEX <- fread("GTEx_gene_median_tpm.gct.gz")
GTEX <- as.data.frame(GTEX)
GTEX <- GTEX[,which(colnames(GTEX) %in% c("Description","Whole Blood"))]

means <- apply(expressed.genes, 1, mean)
means <- na.omit(means)

#looking at arbitry cutoffs
high.means <- means[which(means > 800)]
high <- filter.genes(expressed.genes, names(high.means))
high.GTEX <- GTEX[GTEX$`Whole Blood` > 2000,]

high <- pretty.gene.name(high, as.row.names = TRUE)

# comparing overlap: GSTA4 is interesting 
high.GTEX$Description[which(!high.GTEX$Description %in% high$pretty)]
high$pretty[which(!high$pretty %in% high.GTEX$Description)]

#sorting data first 
GTEX <- filter.out.genes(GTEX, "^MT-",col = 1, by.rownames = FALSE) # removal of mito genes
GTEX.sorted <- sort(GTEX$`Whole Blood`, decreasing = TRUE)
top_hun <- GTEX.sorted[1:100]
top.GTEX <- GTEX[which(GTEX$`Whole Blood` %in% top_hun),]

genes.sorted <- sort(means, decreasing = TRUE)
top_hun_genes <- genes.sorted[1:100]
top.genes <- filter.genes(expressed.genes, names(top_hun_genes))
