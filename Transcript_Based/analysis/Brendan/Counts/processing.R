library(tidyverse)

counts <- fc$counts
counts <- counts[, -which(colnames(counts)%in% c("L6.unmatched.sorted.bam", "L3.unmatched.sorted.bam"))] # remove unmatched
# removes nonunique feature from early hiseq run
colnames(counts)[1:8] <- substring(colnames(counts)[1:8], first = 8) 
colnames(counts)[1:8] <- substring(colnames(counts)[1:8], first = 0, last = 13) 
# removes nonunique feature from all other samples
colnames(counts)[9:ncol(counts)] <- substring(colnames(counts)[9:ncol(counts)], first = 0, last = 9) 

#corplot: samples below 0.85 marked for removal 
cor.plots(counts, heatmap = TRUE, phylo = F)
cor.samples <- cor(counts, method = "spearman")
remove <- names(which(apply(cor.samples, 2, sum) < ncol(cor.samples)*.88))

stat <- fc$stat

#stat data processing
stat <- stat[, -which(colnames(stat)%in% c("L6.unmatched.sorted.bam", "L3.unmatched.sorted.bam"))]# remove unmatched
# removes nonunique feature from early hiseq run
colnames(stat)[2:9] <- substring(colnames(stat)[2:9], first = 8) 
colnames(stat)[2:9] <- substring(colnames(stat)[2:9], first = 0, last = 13) 
# removes nonunique feature from all other samples
colnames(stat)[10:ncol(stat)] <- substring(colnames(stat)[10:ncol(stat)], first = 0, last = 9) 
removal <- stat[which(colnames(stat) %in% c("Status", remove))]

#all samples
med <- stat %>% pivot_longer(col= 2:ncol(stat), 
                                names_to = "sample" , 
                                values_to = "count") %>% 
  filter(count > 0)
#med counts value
med.value <- t %>% group_by(sample) %>% summarize(total = sum(count)) %>%summarize(med = median(total))

#removed samples
t <- removal %>% pivot_longer(col= 2:ncol(removal), 
                              names_to = "sample" , 
                              values_to = "count") %>% 
  filter(count > 0)

             
# ratio of mappings
t %>%  ggplot(aes(sample, count, fill = Status)) + 
  geom_bar(stat = "identity", position = "fill") 

t[nrow(t)+1,] <- med.value$med # addition of med data
t[nrow(t),"sample"] <- "Median"
t <- t %>% mutate(ismed = case_when(sample == "Median" ~ T, sample != "Median" ~ F))
options(scipen=10000) # gets rid of scietific notations

#plotting totals
t %>% group_by(sample,ismed) %>% 
  summarize(total = sum(count)) %>% 
  ggplot(aes(sample, total)) + 
  geom_bar(aes(fill = ismed),
           stat = "identity", 
           show.legend = FALSE) + 
  scale_y_log10() + theme_grey()

#not log scaled
t %>% group_by(sample,ismed) %>% 
  summarize(total = sum(count)) %>% 
  ggplot(aes(sample, total)) + 
  geom_bar(aes(fill = ismed),
           stat = "identity", 
           show.legend = FALSE) + theme_grey()

options(scipen=1)

#sample removal
counts <- counts[,-which(colnames(counts) %in% remove)] # poor cor
pheno.dis <- c("L6.TTAGGC", "L3.TTAGGC")
counts <- counts[,-which(colnames(counts) %in% pheno.dis)] # phenotype discrepancy 

#corplot: samples below 0.85 marked for removal 
cor.plots(counts, heatmap = TRUE, phylo = F)

pheno <- pheno.more.samples
pheno$sample <- gsub("_", ".", pheno$`colnames(expressed.genes)`)
pheno$sample <- gsub("-", ".", pheno$sample)
pheno<-pheno[,-which(colnames(pheno) %in% "colnames(expressed.genes)")]

which(!(pheno$sample %in% colnames(counts))) 
colnames(counts)[which(!(colnames(counts)%in% pheno$sample))]
pheno[nrow(pheno)+1,] <- c("W","III.two","blue",colnames(counts)[which(!(colnames(counts)%in% pheno$sample))])

pheno <- pheno %>% dplyr::select(-pheno.colors)
pheno$sample <- sort(pheno$sample)
counts <- counts[,sort(colnames(counts))]
condition <- factor(pheno$pheno)
rownames(pheno) <- pheno$sample
pheno <- pheno %>% dplyr::select(-sample)

pheno <-pheno %>% mutate(pheno = case_when(pheno == "SS" ~ "S", 
                                   pheno == "S" ~ "S", 
                                   pheno == "W"  ~ "W", 
                                   pheno == "C"~ "C"))




library(edgeR)
library(statmod)

#filtering
cpm_log <- cpm(counts, log = TRUE)
median_log2_cpm <- apply(cpm_log, 1, median)
hist(median_log2_cpm)
expr_cutoff <- -1
abline(v = expr_cutoff, col = "red", lwd = 3)
sum(median_log2_cpm > expr_cutoff)
data_clean <- counts[median_log2_cpm > expr_cutoff, ]

cpm_log <- cpm(data_clean, log = TRUE)
cor.plots(cpm_log,phylo = F)

group <- pheno$pheno
group <- factor(group)
design <- model.matrix(~ 0+group)
colnames(design) <- levels(group)
design






y <- DGEList(counts = counts, group = group)
y$samples

keep <- filterByExpr(y, design)
table(keep)


y <- calcNormFactors(y)
y$sample
y <- y[keep, , keep.lib.sizes=FALSE]

AveLogCPM <- aveLogCPM(y)
hist(AveLogCPM)



pch <- c(0,1,2,15,16,17)
colors <- rep(c("darkgreen", "red", "blue"), 2)
plotMDS(y, col=colors[group], pch=pch[group])
legend("bottomright", legend=levels(group), pch=pch, col=colors, ncol=2)

plotMD(y, column=1)
abline(h=0, col="red", lty=2, lwd=2)

#estimate disp
y <- estimateDisp(y, design, robust = T)

y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)



sqrt(y$common.dispersion) # biological coefficient of variation
plotBCV(y)

#extract top tags
et <- exactTest(y)
results_edgeR <- topTags(et, n = nrow(counts), sort.by = "PValue")
head(results_edgeR$table, n =50)

#smear plot
sum(results_edgeR$table$FDR < .1)
plotSmear(et, de.tags = rownames(results_edgeR)[results_edgeR$table$FDR < .1])
abline(h = c(-2, 2), col = "blue")

fit <- glmFit(y, design, robust=TRUE)
head(fit$coefficients)

#plot for glmQLFit
plotQLDisp(fit)
summary(fit$df.prior)

#c vs s
CvsS <- makeContrasts(S-W, levels=design)
res <- glmLRT(fit, contrast=CvsS)
topTags(res, n = 20)

#c vs w
CvsW <- makeContrasts(C-W, levels=design)
res <- glmLRT(fit, contrast=CvsW, coef = 3)
topTags(res, n = 20)

SvsW <- makeContrasts(S-W, levels=design)
res <- glmLRT(fit, contrast=SvsW)
topTags(res, n = 20)

#s vs c
SvsC <- makeContrasts(S-C, levels=design)
res <- glmLRT(fit, contrast=SvsC)
topTags(res, n = 50)

is.de <- decideTestsDGE(res)
summary(is.de)

plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"),legend="topright")



library(DESeq2)
countTable <- counts

dds <- DESeqDataSetFromMatrix(countData = countTable,
                              colData = pheno,
                              design= ~  pheno)


#rlog transfrom
vs <- vst(dds)

meanSdPlot(assay(vs))

#PCA plot 
pcaData <- DESeq2::plotPCA(vs, intgroup = c("pheno"))
ggplot(pcaData$data)  + 
  geom_text(aes(PC1, PC2, color = pheno), label =pheno$Replicates)+
  ylab(pcaData$labels$y)+
  xlab(pcaData$labels$x)

##########exper


dds <- estimateSizeFactors(dds)
sizeFactors(dds)
#total raw counts
colSums(counts(dds))
## Total number of normalized counts per sample
colSums(counts(dds, normalized=T))
## Plot dispersion estimates
dds <- estimateDispersions(dds)
plotDispEsts(dds)
dds <- DESeq(dds)
res <- results(dds, contrast = c("pheno","S","C"))
res_tableOE <- lfcShrink(dds, coef = 2, res=res, type= "apeglm")
plotMA(res, ylim=c(-5,5))
plotMA(res_tableOE, ylim=c(-2,2))
summary(res_tableOE)

nbinomWaldTest(dds)

results(dds, contrast=c("pheno", "C", "W")) # log fold changes 

resApeT <- lfcShrink(dds, coef=2, type="apeglm", lfcThreshold=1)
plotMA(resApeT, cex=.8)
abline(h=c(-1,1), col="dodgerblue", lwd=2)

###diff
design(dds)
dds <- DESeq(dds)
res <- results( dds )

sig <-res[which(res$pvalue < 0.005),]


idx <- which.min(res$pvalue)
counts(dds)[idx, ]
#Normalization
counts(dds, normalized=TRUE)[ idx, ]

library(biomaRt)
#gene name conversions
#mart =  useDataset("hsapiens_gene_ensembl",mart=ensembl)
ensembl <- useEnsembl(biomart="ensembl", 
                      dataset="hsapiens_gene_ensembl")
genes <- rownames(fc$counts)
#G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
results <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                 filter = "ensembl_gene_id_version", 
                 values = genes, mart = ensembl)
for(i in 1:nrow(results)){
  idx <- grep(results$ensembl_gene_id[i], rownames(counts))
  rownames(counts)[idx] <- results$hgnc_symbol[i]
  
}


rownames(counts) 
