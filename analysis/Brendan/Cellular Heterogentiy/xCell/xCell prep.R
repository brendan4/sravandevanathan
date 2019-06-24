# ONE TIME INSTALL 
devtools::install_github('dviraran/xCell')

#results from R package 
library(xCell)
#must remove unique genes names to match signature file
test <- pretty.gene.name(na.omit(expressed.genes), remove.dups = TRUE, as.row.names = TRUE)
# run x.celll with all genes in sig file 
R.results <- xCellAnalysis(test)

# run the test with fitlered sig file
x.cell.whole.blood <- function(date.set){
  raw.scores = rawEnrichmentAnalysis(as.matrix(data.set),
                                     xCell.data$signatures,
                                     xCell.data$genes)
  
  # cleanup data returned
  colnames(raw.scores) = gsub("\\.1","",colnames(raw.scores))
  raw.scores = aggregate(t(raw.scores)~colnames(raw.scores),FUN=mean)
  rownames(raw.scores) = raw.scores[,1]
  raw.scores = raw.scores[,-1]
  raw.scores = t(raw.scores)
  
  #modify which cells types should be included for accuracy 
  whole.blood <- c("Erythrocytes", "B-cells", "Basophils", "CD4+ memory T-cells", "CD4+ naive T-cells", 
                   "CD4+ T-cells", "CD4+ Tcm", "CD4+ Tem", 'CD8+ T-cells', 'CD8+ Tcm', "CD8+ Tem", "cDC",
                   "Class-switched memory B-cells", "CLP", "CMP","DC", "Eosinophils", "Erythrocytes", "HSC",
                   "iDC", "Macrophages", "Macrophages M1", "Macrophages M2", "Mast cells", 'Memory B-cells',
                   'Monocytes', 'MPP', 'MSC', 'naive B-cells', 'Neutrophils', 'NK cells', 'NKT', 'pDC',"Plasma cells",
                   'Platelets', 'pro B-cells', "Tgd cells", 'Th1 cells', 'Th2 cells', 'Tregs')
  
  cell.types = rownames(raw.scores)
  cell.types.use = whole.blood
  cell.types.use = intersect(rownames(raw.scores),whole.blood)
  
  #transformed scores for final results
  transformed.scores = transformScores(raw.scores[cell.types.use,],xCell.data$spill$fv)
  scores = spillOver(transformed.scores,xCell.data$spill$K)
  #s = y
  A = intersect(colnames(test),colnames(scores))
  return(scores = scores[,A])
}

scores <- x.cell.whole.blood(test)
write.csv(scores, "xCell.csv")

hc <- hclust(as.dist(1-cor(t(na.omit(scores)))))
heatmap(na.omit(scores), Rowv=as.dendrogram(hc))
heatmap(na.omit(scores))

library(ComplexHeatmap)

data("full.pheno.table")
par(mar=c(7,4,4,2)+0.1) 
png(filename='celltype_heatmap.png', width=800, height=750)

# use t() for k means clustering of sample no genes 
Heatmap(scores,
        column_names_side = "bottom",
        row_names_side = "left",
        row_hclust_side = "left",
        row_names_gp=gpar(cex=0.6),
        row_hclust_width = unit(3, "cm"),
        clustering_distance_rows ="euclidean",
        clustering_method_rows = "centroid",
        bottom_annotation = HeatmapAnnotation(as.data.frame(pheno.types),
                                                 which = "column",
                                                 show_legend=TRUE)) # number of clusters you want
graphics.off()


data("pheno")
data("pheno.colors")
pheno.colors <- t(pheno.colors)


#results from web software
setwd("C:\\Users\\brendan\\Documents\\sravandevanathan\\analysis\\Brendan\\Cellular Heterogentiy\\xCell\\Web Outout")
xCell <- read.table("xCell_results.txt", sep = "\t")
xCell <- as.data.frame(xCell)
rownames(xCell) <- xCell[,1]
xCell <- xCell[,-1]
colnames(xCell) <- colnames(expressed.genes)
xCell <- xCell[-1,]


