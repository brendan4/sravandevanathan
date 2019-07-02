library(openxlsx)
setwd("C:/Users/brendan/Documents/sravandevanathan/analysis/Brendan/Cellular Heterogenity/xCell")
xCell.genesig <- read.xlsx("xCell gene signatures.xlsx")


#create pretty names 
test.gene <- rownames(pretty.gene.name(na.omit(expressed.genes), 
                                       as.row.names = TRUE, 
                                       remove.dups = TRUE))

#whole blood cell types
whole.blood <- c("B-cells", "Basophils", "Eosinophils", 
                 "Erythrocytes", 'Neutrophils',"Plasma cells",
                 'Platelets')

#OPTIONAL: subset of gene sig file based on whole blood compostion
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

### for all cell types in the gene signature
filt.xCell.genesig <- xCell.genesig
# checking for lost genes: found in sig not in our pretty data
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

# removing duplicates from the xCell lost table: total # of lot genes
lost.names <- lost[!duplicated(lost$gene),]

# looking for gene in less filtered expressed.genes
found <- filter.genes(expressed.genes, gene.list = lost.names$gene, lazy= FALSE)
# only one found
lost[grep("GPR85", lost$gene),]

#opening GTEX median expression for cell types
library(data.table)
library(ggplot2)
setwd("C:/Users/brendan/Documents/sravandevanathan/analysis/Brendan/Cellular Heterogenity/xCell")
GTEX <- fread("GTEx_gene_median_tpm.gct.gz")
GTEX <- as.data.frame(GTEX)
GTEX <- GTEX[,which(colnames(GTEX) %in% c("Description","Whole Blood"))]
GTEX.sub <- GTEX[which(GTEX$Description %in% lost.names$gene),]

point.labels <- GTEX.sub[which(GTEX.sub$`Whole Blood` > 1), ]

#plotting gene names 
ggplot(data = GTEX.sub, aes(x = GTEX.sub$Description, y = GTEX.sub$`Whole Blood`)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  xlab("Gene")+
  ylab("TPM")+
  geom_point(alpha = .4)+
  geom_text(data = point.labels, aes(x = point.labels$Description, 
                                     y = point.labels$`Whole Blood`, 
                                     label = point.labels$Description))
