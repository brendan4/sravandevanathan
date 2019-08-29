library(propr)

#data and cleaning
setwd("~/sravandevanathan/analysis/Brendan/reads")
reads <- read.delim2("concordant.txt", sep = ",", header = FALSE)
reads <- reads[,c(1,2,3,4)]
colnames(reads) <- c("sample", "concordant", "discordant", "leftovers")
reads <- reads[-which(reads$sample %in% c("L3_unmatc", "L6_unmatc")),] # removal of unmatched

reads <- reads[-which(!reads$sample %in% colnames(expressed.genes.GEN)), ]

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
test <- unnormilize(expressed.genes.GEN, reads)
test[1,1]/reads[which(reads$sample %in% "L2_ACTTGA"),"concordant"] == expressed.genes.GEN[1,1]
test[14970,6]/reads[which(reads$sample %in% colnames(test)[6]),"concordant"] == expressed.genes.GEN[14970,6]

#for clr
test <- na.omit(test)
test <- pretty.gene.name(test, as.row.names = T, remove.dups = T)
test <- t(test)
keep <- apply(test, 2, function(x) sum(x >= 100) >= 10)
phi <- propr(test, metric = "rho", select = keep)

#Subset option 1
phiHBB.HBG1 <- subset(phi, select = DBA.related.Ribo )
look<-phiHBB.HBG1@matrix
