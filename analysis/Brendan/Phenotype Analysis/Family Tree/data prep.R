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

#SKIP: possible mistake in phenotype table indivduals correlate however they do not appear to have same sex 
sub <- expressed.genes[,which(colnames(expressed.genes) %in% c("L6_TTAGGC", "L3_ACTTGA"))]
sex.genes <- c("XIST", "USP9Y", "UTY", "RPS4Y1", "TSIX", "PRKY", 'DDX3Y', 'RPS4Y1')
sub.filt <- filter.genes(sub, sex.genes)
sub.hem <- filter.genes(sub,"hb")
setwd("C:/Users/brendan/Documents/sravandevanathan/analysis/Brendan/Working Dataset/02_Filtered Data")
filt.names <- read.delim("sex_related.tab", sep= "")
sub.filt <- filter.genes(sub, filt.names$geneNames)

cor(na.omit(sub), method = "spearman")

#dropping L6_TTAGGC and L3_TTAGGC
df <- df[-which(df$id %in% c("L6_TTAGGC")),]
df <- df[-which(df$person  == 0), ]


#collapse individuals

df <- df[!duplicated(df$person),]
df$id <- df$person
df <- df[, - which(colnames(df) %in% c("person"))]

#manually entering data

df[which(df$id %in% "I.one"), "dadid"] <- NA
df[which(df$id %in% "I.one"), "momid"] <- NA

# I.one kids
df[which(df$id %in% "II.six"), "dadid"] <- "I.one"
df[which(df$id %in% "II.six"), "momid"] <- NA

df[which(df$id %in% "II.seven"), "dadid"] <- "I.one"
df[which(df$id %in% "II.seven"), "momid"] <- NA

df[which(df$id %in% "II.eight"), "dadid"] <- "I.one"
df[which(df$id %in% "II.eight"), "momid"] <- NA

# entering  II.one to the table ( not in phenotype data)
df[length(df)+1,] <- "II.one"
df$famid <- 1
df[which(df$id %in% "II.one"), "sex"] <- 2

df[which(df$id %in% "II.one"), "dadid"] <- "I.one"
df[which(df$id %in% "II.one"), "momid"] <- NA

#entering third level

#II.one children
df[which(df$id %in% "III.one"), "dadid"] <- NA
df[which(df$id %in% "III.one"), "momid"] <- "II.one"

df[which(df$id %in% "III.two"), "dadid"] <- NA
df[which(df$id %in% "III.two"), "momid"] <- "II.one"

df[which(df$id %in% "III.three"), "dadid"] <- NA
df[which(df$id %in% "III.three"), "momid"] <- "II.one"

df[which(df$id %in% "III.four"), "dadid"] <- NA
df[which(df$id %in% "III.four"), "momid"] <- "II.one"

df[which(df$id %in% "III.five"), "dadid"] <- NA
df[which(df$id %in% "III.five"), "momid"] <- "II.one"

df[which(df$id %in% "III.six"), "dadid"] <- NA
df[which(df$id %in% "III.six"), "momid"] <- "II.one"

#II.8 children

df[which(df$id %in% "III.eleven"), "dadid"] <- NA
df[which(df$id %in% "III.eleven"), "momid"] <- "II.eight"

df[which(df$id %in% "III.twelve"), "dadid"] <- NA
df[which(df$id %in% "III.twelve"), "momid"] <- "II.eight"

df[which(df$id %in% "III.thirteen"), "dadid"] <- NA
df[which(df$id %in% "III.thirteen"), "momid"] <- "II.eight"

df[which(df$id %in% "III.fourteen"), "dadid"] <- NA
df[which(df$id %in% "III.fourteen"), "momid"] <- "II.eight"

# IV.five

df[which(df$id %in% "IV.five"), "dadid"] <- NA
df[which(df$id %in% "IV.five"), "momid"] <- "III.two"


# adding

pedigree(id = df$id, dadid = df$dadid, momid = df$momid, sex = as.numeric(df$sex))
