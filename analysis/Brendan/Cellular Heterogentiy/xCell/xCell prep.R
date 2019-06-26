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

#Heatmaps 

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
        clustering_distance_rows ="maximum",
        clustering_method_rows = "complete",
        bottom_annotation = HeatmapAnnotation(type = full.pheno.table[,c(2)],
                                              col = list(type = c("SS" =  "red", "S" = "yellow", "C"= "grey", "W"="black")),
                                                 which = "column",
                                                 show_legend=TRUE),
        top_annotation = HeatmapAnnotation(type = full.pheno.table[,c(3)],
                                              which = "column",
                                              show_legend=TRUE))
graphics.off()

#annotation setup
data("full.pheno.table")
pheno.table <- full.pheno.table[,c(2,1)]
rownames(pheno.table) <- pheno.table[,2]
pheno.table[,2] <- full.pheno.table[,3]
colnames(pheno.table) <- c("Pheno", "Replicates")

# Specify colors for annnotation 
Pheno = c("red", "yellow", "grey", "black")
names(Pheno) = c("SS", "S", "C", "W")
ann_colors = list(Pheno = Pheno)

cor.plots(scores, method = "spearman", annotation = pheno.table, colors = ann_colors)




#results from web software
setwd("C:\\Users\\brendan\\Documents\\sravandevanathan\\analysis\\Brendan\\Cellular Heterogentiy\\xCell\\Web Outout")
xCell <- read.table("xCell_results.txt", sep = "\t")
xCell <- as.data.frame(xCell)
rownames(xCell) <- xCell[,1]
xCell <- xCell[,-1]
colnames(xCell) <- colnames(expressed.genes)
xCell <- xCell[-1,]


