

#data and cleaning
setwd("~/sravandevanathan/analysis/Brendan/reads")
reads <- read.delim2("concordant.txt", sep = ",", header = FALSE)
reads <- reads[,c(1,2,3,4)]
colnames(reads) <- c("sample", "concordant", "discordant", "leftovers")
reads <- reads[-which(reads$sample %in% c("L3_unmatc", "L6_unmatc")),] # removal of unmatched
reads <- reads[-which(!reads$sample %in% colnames(expressed.genes.GEN)), ]

#unnormalizing the data
unnormilize <- function(data.set, reads){
  for(i in 1:nrow(reads)){
    reads[i,"concordant"]
    match <- which(colnames(data.set) %in% reads[i,"sample"])
    data.set[,match] <- data.set[,match] * reads[i,"concordant"]
  }
  data.set <- data.set/1000000
  return(data.set)
}

# sanity check for the data
unnorm.expressed.genes <- unnormilize(expressed.genes.GEN, reads)
unnorm.expressed.genes[1,1]/
  reads[which(reads$sample %in% "L2_ACTTGA"),"concordant"]*1000000 == expressed.genes.GEN[1,1]
unnorm.expressed.genes[14970,6]/
  reads[which(reads$sample %in% colnames(unnorm.expressed.genes)[6])
        ,"concordant"]*1000000 == expressed.genes.GEN[14970,6]

#for clr
test <- pretty.gene.name(test, as.row.names = T, remove.dups = T)
unnorm.expressed.genes <- t(unnorm.expressed.genes)
keep <- apply(unnorm.expressed.genes, 2, function(x) sum(x >= 100) >= 10)
phi <- propr(unnorm.expressed.genes, metric = "rho", select = keep)

#finding DBA genes
#DBA releated ribo genes
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6416817/
DBA.related.Ribo <- c("RPL5", "RPL11", "RPL35A", "RPS7", "RPS10", "RPS17", "RPS19", "RPS24", "RPS26",
                      "RPL3", "RPL7", "RPL9", "RPL14", "RPL19", "RPL23A", "RPL26", "RPL35", "RPL36", "RPS8"
                      ,"RPS15", "RPS27A", "RPL18")
DBA <- filter.genes(expressed.genes.GEN, DBA.related.Ribo, lazy = F)
DBA <- data.frame(pretty = DBA$pretty, row.names = rownames(DBA)) #simplifing the data

#Subset option 1
phi.DBA <- subset(phi, select = rownames(DBA))
DBA.matrix <- phi.DBA@matrix
snapshot(phi.DBA, prompt = TRUE, plotly = TRUE)

phi.DBA.sub <- subset(phi, select = rownames(DBA)[which(DBA$pretty %in% c("RPL11", "RPS7", "RPS10"))])
smear(phi.DBA.sub, plotly = TRUE)

library(ComplexHeatmap)

# use t() for k means clustering of sample no genes 
Heatmap(DBA.matrix,
        column_names_side = "bottom",
        row_names_side = "left",
        row_hclust_side = "left",
        row_hclust_width = unit(3, "cm"))
graphics.off()