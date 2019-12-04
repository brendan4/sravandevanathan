#comparison of L2_ACAGTG: sample that was thrown out
L2 <- gene_abundance.tab
L2OLD <- gene_abundance

length(L2) == length (L2OLD)

summary(L2$FPKM)
summary(L2OLD$FPKM)
hist(log(L2$FPKM+0.1))
hist(log(L2OLD$FPKM+0.1))

boxplot(log2(L2$FPKM+1) 
        , ylab="log2(FPKM)", 
        main="Distribution of FPKMs for all libraries")

# comparion of L2_CTTGTA
setwd("C:\\Users\\brendan\\Documents\\sravandevanathan\\run_june_3_2019\\ballgown\\L2_CTTGTA")
newL2_CTTGTA <- read.delim("gene_abundance.tab.gz", sep = "\t")
setwd("C:\\Users\\brendan\\Documents\\sravandevanathan\\ballgown\\L2_CTTGTA")
oldL2_CTTGTA <- read.delim("gene_abundance.tab", sep = "\t")

summary(newL2_CTTGTA$FPKM)
summary(oldL2_CTTGTA$FPKM)

hist(log10(newL2_CTTGTA$FPKM+0.1))
hist(log10(oldL2_CTTGTA$FPKM+0.1))



