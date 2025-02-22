library("readxl")
setwd("C:/Users/brendan/Documents/sravandevanathan/")
pheno <- read_excel("061518_phenotypes.xlsx")

# all mutant carriers 
I.one <- c('L3_ACTTGA', "L6_TTAGGC")
II.eight <- c("L3_GGCTAC")
III.two <- c("L3_CGATGT", "L6_ACAGTG", "L2_ATCACG","L2_GCCAAT", "L6_TGACCA")
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
III.twelve <- c("L6_CTTGTA", "L2_TAGCTT")
III.fourteen <- c("L6_ATCACG", "L2_GGCTAC")

carriers <- c(I.one,II.eight,III.two,III.three, III.four)
some.symp <- III.thirteen
severe <- IV.five

pheno <- as.data.frame(colnames(expressed.genes))
pheno$pheno <- 0
pheno[which(pheno[,1] %in% carriers), 2] <- "C"
pheno[which(pheno[,1] %in% some.symp), 2] <- "S"
pheno[which(pheno[,1] %in% severe), 2] <- "SS"
pheno[which(pheno[,2] == 0), 2] <- "W"

#OPTIONAL: colors
pheno$pheno.colors <- 0
pheno[which(pheno[,2] == "C"),3] <- "yellow"
pheno[which(pheno[,2] == "SS"),3] <- "red"
pheno[which(pheno[,2] == "S"),3] <- "orange"
pheno[which(pheno[,2] == "W"),3] <- "blue"
pheno.colors <- t(pheno[,-2])
colnames(pheno.colors) <- pheno.colors[1,]
pheno.colors <- pheno.colors[-1,]

#mutant list 
carriers <- list(I.one = c('L3_ACTTGA', "L6_TTAGGC"), 
                 II.eight = c("L3_GGCTAC"),
                 III.two = c("L3_CGATGT", "L6_ACAGTG", "L2_ATCACG","L2_GCCAAT", "L6_TGACCA"),
                 III.three = c("L6_GCCAAT", "L6_ACTTGA","L2_CAGATC"), 
                 III.four = c("L3_TGACCA", "L6_CAGATC", "L2_CGATGT", "L2_ACTTGA"))

some.symp <- list(III.thirteen = c("L6_GGCTAC", "L2_GATCAG"))
severe <- list(IV.five = c("L2_CTTGTA", "L6_CGATGT"))
mutants <- list(carriers, some.symp, severe)

wildtype <- list(II.six = c("L3_GATCAG"),
                 II.seven = c("L3_TAGCTT"),
                 III.one = c("L3_ATCACG", "L3_CTTGTA", "L2_TGACCA", "L2_ACAGTG"),
                 III.five = c("L3_ACAGTG","L2_TTAGGC"),
                 III.six = c("L3_GCCAAT","L6_GATCAG"),
                 III.eleven = c("L3_CAGATC", "L6_TAGCTT"), 
                 III.twelve = c("L6_CTTGTA", "L2_TAGCTT"), 
                 III.fourteen = c("L6_ATCACG", "L2_GGCTAC"))
full.pheno <- list(mutants, wildtype)

names(full.pheno)<- c("Mutant", "Wildtype")
names(full.pheno$Mutant) <- c("Carrier", "Intermediate", "Severe")


data("pheno")
data("full.pheno")
data("pheno.colors")

# create a table with pheno info and mutation class
full.pheno.table <- pheno 

add.reps <- function(data.set, full.pheno){
  data.set$Replicates <- 0
  data.set$colors <- 0 
  for(person in 1:length(full.pheno$Wildtype)){
    name <- names(full.pheno$Wildtype)[person]
    samples <- full.pheno$Wildtype[[person]]
    data.set[which(pheno[, 1] %in% samples),3] <- name
  }
  for (subcat in 1:length(full.pheno$Mutant)){
    for (person in 1:length(full.pheno$Mutant[[subcat]])){
      name <- names(full.pheno$Mutant[[subcat]])[person]
      samples <- full.pheno$Mutant[[subcat]][[person]]
      data.set[which(pheno[, 1] %in% samples),
            which(colnames(data.set) %in% "Replicates")] <- name
    }
  }
  return(data.set)
}

full.pheno.table <- add.reps(full.pheno.table, full.pheno = full.pheno)

#adding full.pheno.table
full.pheno.table$pheno.colors <- 0
full.pheno.table[which(pheno[,2] == "C"),4] <- "yellow"
full.pheno.table[which(pheno[,2] == "SS"),4] <- "red"
full.pheno.table[which(pheno[,2] == "S"),4] <- "orange"
full.pheno.table[which(pheno[,2] == "W"),4] <- "blue"


