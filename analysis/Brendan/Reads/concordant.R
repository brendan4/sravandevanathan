setwd("~/sravandevanathan/analysis/Brendan/reads")
reads <- read.delim2("concordant.txt", sep = ",", header = FALSE)
concordant <- concordant[,c(1,2)]
colnames(concordant) <- c("smaple", "concordant reads")
