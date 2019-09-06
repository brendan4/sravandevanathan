library(ComplexHeatmap)

hemo = c("HBA", "HBB", "HBD", "HBG","HBZ")
hemo = c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")

hemo.genes <- filter.genes(expressed.genes.GEN, hemo)
hemo.genes <- pretty.gene.name(hemo.genes, as.row.names = TRUE)

data("pheno.more.samples")
par(mar=c(7,4,4,2)+0.1) 
png(filename='top_heatmap.png', width=1000, height= 800)


Heatmap(log(hemo.genes+0.1),
        column_names_side = "bottom",
        column_title_gp = gpar(cex=.50),
        row_names_side = "left",
        row_hclust_side = "left",
        row_names_gp=gpar(cex=0.75),
        row_hclust_width = unit(3, "cm"),
        clustering_distance_columns = "spearman",
        clustering_method_columns = "complete",
        clustering_distance_rows ="maximum",
        clustering_method_rows = "complete",
        bottom_annotation = HeatmapAnnotation(type = pheno.more.samples[,c(2)],
                                              col = list(type = c("SS" =  "red", "S" = "yellow", "C"= "grey", "W"="black")),
                                              which = "column",
                                              show_legend=TRUE))
graphics.off()


