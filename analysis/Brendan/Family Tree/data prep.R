# adding person data
library(kinship2)
full.pheno.table$Replicates
full.pheno.table <- data.frame(full.pheno.table, stringsAsFactors = FALSE)

df <- data.frame(id = colnames(expressed.genes), sex = 0, dadid = 0, momid = 0, famid = 1, stringsAsFactors = FALSE)
df$person = 0
for(row in 1:nrow(full.pheno.table)){
  ids <-  full.pheno.table[row,3]
  samples <- full.pheno.table[which(full.pheno.table[,3] %in% ids),1]
  df[which(df$id %in% samples), 6] <- ids
}

#adding sex data
data("sex")
sex$Person <- rownames(sex)
names <- rownames(sex[which(sex$Gender %in% "M"),])
df[which(df$id %in% names), "sex"] <- 1

names <- rownames(sex[which(sex$Gender %in% "F"),])
df[which(df$id %in% names), "sex"] <- 2

sub <- expressed.genes[,which(colnames(expressed.genes) %in% c("L6_TTAGGC", "L3_ACTTGA"))]
sub.xist <- sub[grep("^XIST", rownames(sub)), ]

