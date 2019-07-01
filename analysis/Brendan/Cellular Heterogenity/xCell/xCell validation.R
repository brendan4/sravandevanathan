
library(openxlsx)
setwd("C:/Users/brendan/Documents/sravandevanathan/analysis/Brendan/Cellular Heterogenity/xCell")
xCell.genesig <- read.xlsx("xCell gene signatures.xlsx")



test.gene <- rownames(pretty.gene.name(na.omit(expressed.genes), 
                                       as.row.names = TRUE, 
                                       remove.dups = TRUE))
lost <- data.frame(0)
for(row in 1:nrow(xCell.genesig)){
  print(paste("Checking celltype:", xCell.genesig[row, 1]))
  
  for(gene in 1:xCell.genesig[row, 2]){
    
    if (xCell.genesig[row, gene +2] %in% test.gene){
      print(paste(xCell.genesig[row, gene +2],": found in expressed.genes"))
    }else{
      print(paste(xCell.genesig[row, gene +2],": NOT found in expressed.genes"))
    }
  }
}
