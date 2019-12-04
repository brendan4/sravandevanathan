
setwd("~/sravandevanathan/counts_script")
counts <- read.table("unprocessed_counts.tab", row.names = 1, check.names = F)

# clean sample names
colnames(counts) <- gsub(".sorted.bam$", "",
                         colnames(counts)) # removes files extension
colnames(counts) <- gsub("^121317.", "",
                         colnames(counts)) # removes tag from lib samples

# removal of samples
cor.table <- cor(counts, method = "spearman")
pheno.dis <- c("L6.TTAGGC", "L3.TTAGGC")
remove <- c(names(which(apply(cor.table, 2, sum) < ncol(cor.samples)*.88)), pheno.dis)
counts <- counts %>% 
  dplyr::select(-contains("unmatched")) %>% 
  dplyr::select(-one_of(remove))

setwd("~/sravandevanathan/counts_script/Output/Data")
write.table(counts, "counts.tab", sep = "\t")

