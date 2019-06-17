library("readxl")
setwd("C:/Users/brendan/Documents/sravandevanathan/")
pheno <- read_excel("061518_phenotypes.xlsx")

carriers <- c(I.one,II.eight,III.two,III.three, III.four)
some.symp <- III.thirteen
severe <- IV.five

pheno <- as.data.frame(colnames(expressed.genes))
pheno$pheno <- 0
pheno[which(pheno[,1] %in% carriers), 2] <- "C"
pheno[which(pheno[,1] %in% some.symp), 2] <- "S"
pheno[which(pheno[,1] %in% severe), 2] <- "SS"
pheno[which(pheno[,2] == 0), 2] <- "W"

# all mutant carriers 
I.one <- c('L3_ACTTGA', "L6_TTAGGC")
II.eight <- c("L3_GGCTAC")
III.two <- c("L3_CGATGT", "L6_ACAGTC", "L2_ATCACG","L2_GCCAAT", "L6_TGACCA")
III.three <- c("L6_GCCAAT", "L6_ACTTGA","L2_CAGATC")
III.four <- c("L3_TGACCA", "L6_CAGATC", "L2_CGATGT", "L2_ACTTGA")
III.thirteen <- c("L6_GGCTAC", "L2_GATCAG")
IV.five <- c("L2_CTTGTA", "L6_CGATGT")

# wildtype 
II.six <- c("L3_GATCAG")
II.seven <- c("L3_TAGCTT")
III.one <- c("L3_ATCACG", "L3_CTTGTA", "L2_TGACCA", "L2_ACAGTG")
III.five <- c("L3_ACAGTG","L2_TTAGGC")
III.six <- c("L3_GCCAAT","L6_GATCAG")
III.eleven <- c("L3_CAGATC", "L6_TAGCTT")
III.twelve <- C("L6_CTTGTA", "L2_TAGCTT")
III.fourteen <- c("L6_ATCACG", "L2_GGCTAC")


