library(tidyverse)
# importing data
setwd("~/sravandevanathan/counts_script")
counts <- read.table("Input_data/unprocessed_counts.tab", row.names = 1, check.names = F)
pheno <- read.table("Input_data/pheno.tab",row.names = 1, check.names = F)

# clean sample names
colnames(counts) <- gsub(".sorted.bam$", "",
                         colnames(counts)) # removes files extension
colnames(counts) <- gsub("^121317.", "",
                         colnames(counts)) # removes tag from lib samples

# removal of samples
counts <- counts %>% dplyr::select(-contains("unmatched"))
cor.table <- cor(counts, method = "spearman")
pheno.dis <- c("L6.TTAGGC", "L3.TTAGGC")
remove <- c(names(which(apply(cor.table, 2, sum) < ncol(cor.samples)*.88)), pheno.dis)
counts <- counts %>% dplyr::select(-one_of(remove))

#WARN: saving options
setwd("~/sravandevanathan/counts_script/Output/Data")
write.table(counts, "counts.tab", sep = "\t")

# gene name conversions
library(biomaRt)
genes <- rownames(counts)
genes <- gsub("\\.\\d*$", "", genes)
ensembl <- useEnsembl(biomart="ensembl", 
                      dataset="hsapiens_gene_ensembl")

results <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                 filter = "ensembl_gene_id", 
                 values = genes, mart = ensembl)
# replacing rownames in counts matrix
for(i in 1:nrow(results)){
  idx <- grep(paste(results$ensembl_gene_id[i], "\\.\\d*$", sep= ""), rownames(counts)) # finds match
  if (nchar(results$hgnc_symbol[i]) == 0){ # no gene id replaced with ensembl id
    rownames(counts)[idx] <- results$ensembl_gene_id[i]}
  else if (results$hgnc_symbol[i] %in% rownames(counts)){ # prevents duplicates gene names
    rownames(counts)[idx] <- results$ensembl_gene_id[i]
  }else{ # replacing with gen id
    rownames(counts)[idx] <- results$hgnc_symbol[i]
  }
}

#DESeq2 analysis
library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = counts.t,
                              colData = pheno,
                              design = ~  pheno) # experiment declaration 

# PCA plot with VST transform 
vsd <- vst(dds)
plotPCA(vsd, intgroup = c("pheno"))

pcaData <- DESeq2::plotPCA(vsd, intgroup = c("pheno"))
ggplot(pcaData$data)  + # with replicate labels
  geom_text(aes(PC1, PC2, color = pheno), label = pheno$Replicates)+
  ylab(pcaData$labels$y)+
  xlab(pcaData$labels$x)
ggplot(pcaData$data)  +  # without replicate labels 
  geom_text(aes(PC1, PC2, color = pheno), label = rownames(pheno))+
  ylab(pcaData$labels$y)+
  xlab(pcaData$labels$x)

# differential expression 
dds <- DESeq(dds)
plotDispEsts(dds)
res <- results(dds, contrast = c("pheno", "C", "W")) #can specify group comparisons
summary(res)

# top results
top <-rownames(res)[which(res$padj < 0.00005)]
plotCounts(dds, "ENSG00000142676.14", "pheno") # plot normalized counts of specific genes by group
plotMA(res, ylim=c(-4,4))

# skrinkage of effect size: mainly for plotting
resultsNames(dds) # select coef for srinkage
resLFC <- lfcShrink(dds, coef=4, type="apeglm")
plotMA(resLFC, ylim=c(-2,2))

# saving names of the significant genes in MA plot
idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]
