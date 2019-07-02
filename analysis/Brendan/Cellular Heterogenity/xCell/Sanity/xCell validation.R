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
gene.sig.filt <- function(gene.sig, cell.types){
  corr.counter <- 0 
  filt.genesig <- data.frame()
  
  for(i in 1:length(cell.types)){
    pattern = grep(paste("^",cell.types[i], sep= ""), gene.sig$Celltype_Source_ID)
    
    if (length(pattern) == 0 ){
      print(paste(cell.types[i]), ": not found")
    } else {
      corr.counter <- corr.counter + 1
    }
    filt.genesig <- rbind(filt.genesig, gene.sig[pattern, ]) 
  }
  return(filt.genesig)
}

filt.xCellgenesig <- gene.sig.filt(gene.sig = xCell.genesig, cell.types = whole.blood)

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

setwd("C:/Users/brendan/Documents/sravandevanathan/analysis/Brendan/Cellular Heterogenity/xCell/Sanity/Datasets")
GTEX <- fread("GTEx_gene_median_tpm.gct.gz")
GTEX <- as.data.frame(GTEX)
GTEX <- GTEX[,which(colnames(GTEX) %in% c("Description","Whole Blood"))]

GTEX.sub <- GTEX[which(GTEX$Description %in% lost.names$gene),]
point.labels <- GTEX.sub[which(GTEX.sub$`Whole Blood` > .1), ]

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

# filtering all genes from wholeblood subset 
sig.all.genes <- data.frame(gene = 0 , celltype = 0)
gene.counter <- 0 
filt.genesig <- gene.sig.filt(gene.sig = xCell.genesig, cell.types = whole.blood)

for(row in 1:nrow(filt.genesig)){
  num.genes <- filt.genesig[row,2]
  genes <- filt.genesig[row, 3:(num.genes+2)]

  for(gene in 1:length(genes)){
    gene.counter <- gene.counter + 1 
    sig.all.genes[gene.counter,"gene"] <- filt.genesig[row, gene +2]
    sig.all.genes[gene.counter, "celltype"] <- filt.genesig[row, 1]
  }
}

#duplicated handling
sum(duplicated(sig.all.genes$gene))
sig.all.genes <- sig.all.genes[-which(duplicated(sig.all.genes$gene)),]

# are any genes in sig.all not found in GTEX data?
sum(!which(sig.all.genes$gene %in% GTEX$Description))

# plotting og GTEX data subseted my sig.all.genes
GTEX.sub <- GTEX[which(GTEX$Description %in% sig.all.genes$gene),]
GTEX.sub$`Whole Blood` <- log(GTEX.sub$`Whole Blood`+0.001)
point.labels <- GTEX.sub[which(GTEX.sub$`Whole Blood` > log(0.1)), ]

#plotting gene names 
ggplot(data = GTEX.sub, aes(x = GTEX.sub$Description, y = GTEX.sub$`Whole Blood`)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  xlab("Gene")+
  ylab("ln(TPM)")+
  geom_point(alpha = .4)+
  geom_point(data = point.labels, shape=23, fill="red",
             aes(x = point.labels$Description, 
                                     y = point.labels$`Whole Blood`))+ 
  annotate("text", x = 500, y = 8, label= paste(nrow(point.labels),": above 0.1 TPM")) + 
  annotate("text", x = 500, y= -5.5, label = paste((nrow(GTEX.sub) - nrow(point.labels)),": below 0.1 TPM"))


plot <- ggplot(data = GTEX.sub, aes(x = GTEX.sub$Description, y = GTEX.sub$`Whole Blood`)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  xlab("Gene")+
  ylab("ln(TPM)")+
  geom_point(alpha = .4)

colors <- c("red", "blue", "Yellow", "greenyellow", "green", "darkred", "dodgerblue")
counter <- 0 
cell.list <- list()
for(cell.type in whole.blood){
  counter <- counter + 1
  cell.specific <- sig.all.genes[grep(cell.type, sig.all.genes$celltype),]
  print(cell.type)
  cell.specific <- GTEX.sub[which(GTEX.sub$Description %in% cell.specific$gene), ]
  point.labels <- cell.specific[which(cell.specific$`Whole Blood` > log(0.1)), ]
  assign(cell.type, point.labels)
  print(B-cells)
  print(paste(cell.type,".data", sep = ""))
  cell.list[[counter]] <- cell.type
  plot <- plot + geom_point(data = cell.type, shape=23, fill = colors[counter],
                            aes(x = cell.type$Description, 
                                y = cell.type$`Whole Blood`))
}

plot

# ploting all genes all cell types 
filt.genesig <- xCell.genesig

for(row in 1:nrow(filt.genesig)){
  num.genes <- filt.genesig[row,2]
  genes <- filt.genesig[row, 3:(num.genes+2)]
  
  for(gene in 1:length(genes)){
    gene.counter <- gene.counter + 1 
    sig.all.genes[gene.counter,"gene"] <- filt.genesig[row, gene +2]
    sig.all.genes[gene.counter, "celltype"] <- filt.genesig[row, 1]
  }
}

#duplicated handling
sum(duplicated(sig.all.genes$gene))
sig.all.genes <- sig.all.genes[-which(duplicated(sig.all.genes$gene)),]

# are any genes in sig.all not found in GTEX data?
sum(!which(sig.all.genes$gene %in% GTEX$Description))

# plotting og GTEX data subseted my sig.all.genes
GTEX.sub <- GTEX[which(GTEX$Description %in% sig.all.genes$gene),]
GTEX.sub$`Whole Blood` <- log(GTEX.sub$`Whole Blood`+0.001)
point.labels <- GTEX.sub[which(GTEX.sub$`Whole Blood` > log(0.1)), ]

#plotting gene names

ggplot(data = GTEX.sub, aes(x = GTEX.sub$Description, y = GTEX.sub$`Whole Blood`)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  xlab("Gene")+
  ylab("ln(TPM)")+
  geom_point(alpha = .4)+
  geom_point(data = point.labels, shape=23, fill="red",
             aes(x = point.labels$Description, 
                 y = point.labels$`Whole Blood`))+ 
  annotate("text", x = 2500, y = 8, label= paste(nrow(point.labels),": above 0.1 TPM")) + 
  annotate("text", x = 2500, y= -5.5, label = paste((nrow(GTEX.sub) - nrow(point.labels)),": below 0.1 TPM"))

