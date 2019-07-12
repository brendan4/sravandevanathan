#opening GTEX median expression for cell types: data loadings
library(data.table)
library(ggplot2)

#data generation
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

# comparing nonoverlap: GSTA4 is interesting 
high.GTEX$Description[which(!high.GTEX$Description %in% high$pretty)]
high$pretty[which(!high$pretty %in% high.GTEX$Description)]

#sorting data
#GTEX
GTEX <- filter.out.genes(GTEX, "^MT-",col = 1, by.rownames = FALSE) # removal of mito genes
GTEX.sorted <- sort(GTEX$`Whole Blood`, decreasing = TRUE)
top_hun <- GTEX.sorted[1:100]
top.GTEX <- GTEX[which(GTEX$`Whole Blood` %in% top_hun),]

#expresed.genes means 
genes.sorted <- sort(means, decreasing = TRUE)
top_hun_genes <- genes.sorted[1:100]
top.genes <- filter.genes(expressed.genes, names(top_hun_genes))
top.genes <- pretty.gene.name(top.genes, as.row.names = TRUE)

#nonintersections between the two 
top.GTEX$Description[which(!top.GTEX$Description %in% top.genes$pretty)]
top.genes$pretty[which(!top.genes$pretty %in% top.GTEX$Description)]

# look into genes in GTEX top not in top.genes 
test <- filter.genes(expressed.genes, gene.list = top.GTEX$Description[which(!top.GTEX$Description %in% top.genes$pretty)], lazy = FALSE)

# look into genes in top.genes not in GTEX
test <- filter.genes(expressed.genes, top.genes$pretty[which(!top.genes$pretty %in% top.GTEX$Description)], lazy = FALSE)
test <- filter.out.genes(test, "-")

#intersections between the two
top.GTEX$Description[which(top.GTEX$Description %in% top.genes$pretty)]

#interesting non-intersecting protien 
expressed.genes[grep("^HBD", rownames(expressed.genes)),]


#heatmaps for expressed genes top expresssion
genes.sorted <- sort(means, decreasing = TRUE)
top_hun_genes <- genes.sorted[1:50]
top.genes <- filter.genes(expressed.genes, names(top_hun_genes))
top.genes <- pretty.gene.name(top.genes, as.row.names = TRUE, remove.dups = TRUE)

library(ComplexHeatmap)

data("full.pheno.table")
par(mar=c(7,4,4,2)+0.1) 
png(filename='top_heatmap.png', width=1000, height= 800)


Heatmap(log(top.genes+0.1),
        column_names_side = "bottom",
        row_names_side = "left",
        row_hclust_side = "left",
        row_names_gp=gpar(cex=0.6),
        row_hclust_width = unit(3, "cm"),
        clustering_distance_columns = "spearman",
        clustering_method_columns = "complete",
        clustering_distance_rows ="maximum",
        clustering_method_rows = "complete",
        bottom_annotation = HeatmapAnnotation(type = full.pheno.table[,c(2)],
                                              col = list(type = c("SS" =  "red", "S" = "yellow", "C"= "grey", "W"="black")),
                                              which = "column",
                                              show_legend=TRUE))
graphics.off()

#intersecting genes 
intersection <- top.GTEX$Description[which(top.GTEX$Description %in% top.genes$pretty)]
intersect <- filter.genes(expressed.genes, intersection, lazy = FALSE)
rownames(intersect) <- intersect$pretty
intersect <- intersect[,-length(intersect)]


Heatmap(log(intersect),
        column_names_side = "bottom",
        row_names_side = "left",
        row_hclust_side = "left",
        row_names_gp=gpar(cex=0.6),
        row_hclust_width = unit(3, "cm"),
        clustering_distance_columns = "spearman",
        clustering_method_columns = "complete",
        clustering_distance_rows ="maximum",
        clustering_method_rows = "complete",
        bottom_annotation = HeatmapAnnotation(type = full.pheno.table[,c(2)],
                                              col = list(type = c("SS" =  "red", "S" = "yellow", "C"= "grey", "W"="black")),
                                              which = "column",
                                              show_legend=TRUE))
graphics.off()
