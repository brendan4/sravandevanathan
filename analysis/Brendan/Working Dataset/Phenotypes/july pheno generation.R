data("pheno")
data("full.pheno.table")
july.pheno <- data.frame("one" = 0, "two" = 0)
colnames(july.pheno) <- colnames(pheno)

for (i in 1:length(july)){
  july.pheno[i, 1] <- colnames(july)[i]
}
july.pheno[,3:4] <- 0
colnames(july.pheno) <- colnames(full.pheno.table)

july.pheno[grep("LIB1", july.pheno$`colnames(expressed.genes)`),2] <- "C"
july.pheno[grep("LIB1", july.pheno$`colnames(expressed.genes)`),3] <- "I.one"
july.pheno[grep("LIB1", july.pheno$`colnames(expressed.genes)`),4] <- "yellow"

july.pheno[grep("LIB3", july.pheno$`colnames(expressed.genes)`),2] <- "W"
july.pheno[grep("LIB3", july.pheno$`colnames(expressed.genes)`),3] <- "II.seven"
july.pheno[grep("LIB3", july.pheno$`colnames(expressed.genes)`),4] <- "blue"

july.pheno[grep("LIB5", july.pheno$`colnames(expressed.genes)`),2] <- "W"
july.pheno[grep("LIB5", july.pheno$`colnames(expressed.genes)`),3] <- "III.one"
july.pheno[grep("LIB5", july.pheno$`colnames(expressed.genes)`),4] <- "blue"

july.pheno[grep("LIB6", july.pheno$`colnames(expressed.genes)`),2] <- "C"
july.pheno[grep("LIB6", july.pheno$`colnames(expressed.genes)`),3] <- "II.two"
july.pheno[grep("LIB6", july.pheno$`colnames(expressed.genes)`),4] <- "yellow"


all.pheno.data <- rbind(full.pheno.table, july.pheno)
