library("readxl")
setwd("C:/Users/brendan/Documents/sravandevanathan/")
pheno <- read_excel("061518_phenotypes.xlsx")

carriers <- c("L3_CGATGT", "L6_ACAGTC", "L2_ATCACG","L2_GCCAAT","L6_GCCAAT",
              "L6_CAGATC","L2_TGACCA","L2_ACAGTG","L2_GCCAAT","L2_CAGATC","L2_ACTTGA")
some.symp <- c("L6_GGCTAC", "L2_GATCAG")
severe <- c("L6_CGATGT", "L2_CTTGTA")

pheno <- as.data.frame(colnames(expressed.genes))
pheno$pheno <- 0
pheno[which(pheno[,1] %in% carriers), 2] <- "C"
pheno[which(pheno[,1] %in% some.symp), 2] <- "S"
pheno[which(pheno[,1] %in% severe), 2] <- "SS"
pheno[which(pheno[,2] == 0), 2] <- "W"

# all mutant carriers 
I.one <- c('L3_ACTTGA', "L6_TTAGGC")
II.eigth <- c("L3_GGCTAC")
III.two <- c("L3_CGATGT", "L6_ACAGTC", "L2_ATCACG","L2_GCCAAT", "L6_TGACCA")
III.three <- c("L6_GCCAAT", "L6_ACTTGA","L2_CAGATC")
III.four <- c("L3_TGACCA", "L6_CAGATC", "L2_CGATGT", "L2_ACTTGA")
III.thirteen <- c("L6_GGCTAC", "L2_GATCAG")
IV.five <- c("L2_CTTGTA", "L6_CGATGT")






