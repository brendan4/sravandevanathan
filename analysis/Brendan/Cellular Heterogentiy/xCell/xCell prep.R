# ONE TIME INSTALL 
devtools::install_github('dviraran/xCell')

#results from R package 
library(xCell)
test <- pretty.gene.name(na.omit(expressed.genes), remove.dups = TRUE, as.row.names = TRUE)
R.results <- xCellAnalysis(test)

  # run the test with 64 sig file 
  raw.scores = rawEnrichmentAnalysis(as.matrix(test),
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
  scores = scores[,A]
  
write.csv(scores, "xCell.csv")
  

#results from web software
setwd("C:\\Users\\brendan\\Documents\\sravandevanathan\\analysis\\Brendan\\Cellular Heterogentiy\\xCell\\Web Outout")
xCell <- read.table("xCell_results.txt", sep = "\t")
xCell <- as.data.frame(xCell)
rownames(xCell) <- xCell[,1]
xCell <- xCell[,-1]
colnames(xCell) <- colnames(expressed.genes)
xCell <- xCell[-1,]
