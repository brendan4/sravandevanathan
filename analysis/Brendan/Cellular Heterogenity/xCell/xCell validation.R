
library(openxlsx)
setwd("C:/Users/brendan/Documents/sravandevanathan/analysis/Brendan/Cellular Heterogenity/xCell")
xCell.genesig <- read.xlsx("xCell gene signatures.xlsx")



test.gene <- rownames(pretty.gene.name(na.omit(expressed.genes), 
                                       as.row.names = TRUE, 
                                       remove.dups = TRUE))

whole.blood <- c("B-cells", "Basophils", "CLP", "CMP","DC", "Eosinophils", 
                 "Erythrocytes", "Macrophages", 'Neutrophils',"Plasma cells",
                 'Platelets')
xCell.genesig[which(xCell.genesig$Celltype_Source_ID %in% whole.blood), ]

for(i in 1:length(whole.blood)){
  xCell.genesig[grep(paste("^",whole.blood[i], sep= ""), xCell.genesig$Celltype_Source_ID)]
}
lost <- c()
for(row in 1:nrow(xCell.genesig)){
  print(paste("Checking celltype:", xCell.genesig[row, 1]))
  
  for(gene in 1:xCell.genesig[row, 2]){
    
    if (xCell.genesig[row, gene +2] %in% test.gene){
      print(paste(xCell.genesig[row, gene +2],": found in expressed.genes"))
    }else{
      print(paste(xCell.genesig[row, gene +2],": NOT found in expressed.genes"))
      lost <- c(lost, xCell.genesig[row, gene +2])
    }
  }
}


view <- filter.genes(expressed.genes,gene.list = lost$lost, lazy= FALSE)
