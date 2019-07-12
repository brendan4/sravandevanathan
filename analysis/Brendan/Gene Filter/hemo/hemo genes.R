hemo = c("HBA", "HBB", "HBD", "HBG","HBZ")
hemo = c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")

hemo.genes <- filter.genes(expressed.genes, hemo)

data("full.pheno.table")
par(mar=c(7,4,4,2)+0.1) 
png(filename='top_heatmap.png', width=1000, height= 800)


Heatmap(log(hemo.genes+0.1),
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