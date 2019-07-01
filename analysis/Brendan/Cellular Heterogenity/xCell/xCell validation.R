
library(openxlsx)
setwd("C:/Users/brendan/Documents/sravandevanathan/analysis/Brendan/Cellular Heterogenity/xCell")
xCell.genesig <- read.xlsx("xCell gene signatures.xlsx")



test.gene <- rownames(pretty.gene.name(na.omit(expressed.genes), 
                                       as.row.names = TRUE, 
                                       remove.dups = TRUE))

whole.blood <- c("B-cells", "Basophils", "CLP", "CMP","DC", "Eosinophils", 
                 "Erythrocytes", "Macrophages", 'Neutrophils',"Plasma cells",
                 'Platelets')

#fitler of gene sig file based on whole blood compostion
filt.xCell.genesig <- data.frame()
corr.counter = 0 
for(i in 1:length(whole.blood)){
  pattern = grep(paste("^",whole.blood[i], sep= ""), xCell.genesig$Celltype_Source_ID)
  if (length(pattern) == 0 ){
    print(paste(whole.blood[i]), ": not found")
  } else {
    corr.counter <- corr.counter + 1
  }
  filt.xCell.genesig <- rbind(filt.xCell.genesig, xCell.genesig[pattern, ]) 
}

# checking for lost genes: found in sig not in our data
lost <- data.frame(gene = 0 , celltype = 0)
lost.counter <- 0
for(row in 1:nrow(filt.xCell.genesig)){
  print(paste("Checking celltype:", filt.xCell.genesig[row, 1]))
  
  for(gene in 1:filt.xCell.genesig[row, 2]){
    
    if (filt.xCell.genesig[row, gene +2] %in% test.gene){
      print(paste(filt.xCell.genesig[row, gene +2],": found in expressed.genes"))
    }else{
      lost.counter <- lost.counter + 1
      print(paste(filt.xCell.genesig[row, gene +2],": NOT found in expressed.genes"))
      lost[lost.counter,"gene"] <- filt.xCell.genesig[row, gene +2]
      lost[lost.counter, "celltype"] <- filt.xCell.genesig[row, 1]
    }
  }
}


view <- filter.genes(expressed.genes,gene.list = lost$lost, lazy= FALSE)
expressed.genes[grep("^P",expressed.genes)]

filter.genes(expressed.genes, "KPYR")
