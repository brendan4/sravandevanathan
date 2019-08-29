library(ggplot2)

#data and cleaning
setwd("~/sravandevanathan/analysis/Brendan/reads")
reads <- read.delim2("concordant.txt", sep = ",", header = FALSE)
reads <- reads[,c(1,2,3,4)]
colnames(reads) <- c("sample", "concordant", "discordant", "leftovers")
reads <- reads[-which(reads$sample %in% c("L3_unmatc", "L6_unmatc")),] # removal of unmatched


ggplot(data = reads, aes(x = sample, y = concordant)) + 
  geom_bar(stat = "identity")

concordant <- reads[,c(1,2)]
ggplot(data = concordant, aes( y = concordant)) + 
  geom_boxplot()

removed <- reads[which(!(reads$sample %in% colnames(expressed.genes.GEN))),]
kept <- reads[which((reads$sample %in% colnames(expressed.genes.GEN))),]
kept.means <- apply(kept[,2:4], 2, mean)
removed$sample <- levels(droplevels(removed$sample))
means.row <- nrow(removed)+1
removed[means.row, "sample"] <- "mean"
removed[means.row, 2:4] <- kept.means
removed[,"status"] <- "removed"
removed[means.row, "status"] <- "kept"

ggplot(data = removed, aes(x = sample, y = concordant, color = status)) + 
  geom_bar(stat = "identity", fill = "white")+ 
  theme_minimal()

reads.prep <- reads 
reads.prep[,"noncordant"] <- reads.prep[3] + reads.prep[4] + reads.prep[2]
reads.prep[,"ratio"] <- reads.prep[2]/reads.prep["noncordant"]
