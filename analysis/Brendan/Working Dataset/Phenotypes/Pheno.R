carriers <- c("L3_ATCACG", "L3_CGATGT", "L3_TGACCA","L6_ACAGTG","L6_GCCAAT",
              "L6_CAGATC","L2_TGACCA","L2_ACAGTG","L2_GCCAAT","L2_CAGATC","L2_ACTTGA")
some.symp <- c("L6_GGCTAC", "L2_GATCAG")
severe <- c("L6_CGATGT", "L2_CTTGTA")

pheno <- as.data.frame(colnames(expressed.genes))
pheno$pheno <- 0
pheno[which(pheno[,1] %in% carriers), 2] <- "C"
pheno[which(pheno[,1] %in% some.symp), 2] <- "S"
pheno[which(pheno[,1] %in% severe), 2] <- "SS"
pheno[which(pheno[,2] == 0), 2] <- "W"








