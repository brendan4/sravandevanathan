library(Rsubread)
library(biomaRt)
main <- "/media/brendan/Elements/dba/DBA_121317"
setwd(main)
all.folders <- list.files()
bam.files <- c()
for(i in 1:length(all.folders)){
  setwd(paste(main, all.folders[i],sep = "/"))
  bam.files <- c(bam.files, paste(main, all.folders[i],list.files(pattern = ".bam$"), sep = "/"))
  setwd(main)
}
setwd("..")
main2 <-"/media/brendan/Elements/dba/main_batch"
setwd(main2)
all.folders <- list.files()
for(i in 1:length(all.folders)){
  setwd(paste(main2, all.folders[i],sep = "/"))
  bam.files <- c(bam.files, paste(main2, all.folders[i],list.files(pattern = ".bam$"), sep = "/"))
  setwd(main2)
}

#smallscale testing
#setwd(paste(main2, all.folders[10],sep = "/"))
#bam.files <- c(bam.files, paste(main2, all.folders[10],list.files(pattern = ".bam$"), sep = "/"))

setwd("..")
setwd("annotation")
anno <- list.files(pattern = ".gz$", full.names = TRUE)
fc <- featureCounts(bam.files, annot.ext= anno, isGTFAnnotationFile = T,isPairedEnd=TRUE)
  look <- fc$annotation
  
#gene name conversions
#mart =  useDataset("hsapiens_gene_ensembl",mart=ensembl)
ensembl <- useEnsembl(biomart="ensembl", 
                      dataset="hsapiens_gene_ensembl")
genes <- rownames(fc$counts)
#G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
results <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                 filter = "ensembl_gene_id_version", 
                 values = genes, mart = ensembl)

filters <- listFilters(ensembl)
filters
stat <- fc$counts
(fc$counts)
write.csv(stat, "counts.csv")
write.table(x=data.frame(fc$annotation[,c("GeneID", "Length")], fc$counts, stringsAsFactor = F), file = "counts.txt", quote= F, sep = "\t", row.names = F)

