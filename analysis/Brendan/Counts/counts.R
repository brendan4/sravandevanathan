indir <- "C:/Users/brendan/Documents/RPL11_BAM_SAM/BAM"
list.files(indir)
data("pheno.more.samples")

filenames <- file.path(indir, paste0(pheno.more.samples$`colnames(expressed.genes)`, "_RPL11.bam"))
file.exists(filenames)

library("Rsamtools")
bamfiles <- BamFileList(filenames, yieldSize=2000000)
seqinfo(bamfiles[1])

library("GenomicFeatures")
gtffile <- file.path(indir,"Gencode.gtf")
txdb <- makeTxDbFromGFF(gtffile, format = "gtf", circ_seqs = character())
txdb

ebg <- exonsBy(txdb, by="gene")
ebg

library("GenomicAlignments")
library("BiocParallel")
register(SerialParam())

se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )

#exploring se 
dim(se)
assayNames(se)
head(assay(se), 3)
colSums(assay(se))
rowRanges(se)
str(metadata(rowRanges(se)))
colData(se)
colData(se) <- DataFrame(pheno.more.samples)
colData(se)

se$colnames.expressed.genes.
round( colSums(assay(se)) / 1e6, 1 ) # millions of fragments that uniquely align

library("DESeq2")
dds <- DESeqDataSet(se, design = ~ pheno)

countdata <- assay(se)
head(countdata, 3)

coldata <- colData(se)
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = ~ pheno)


nrows(dds)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
nrow(dds)

# at least 3 samples with a count of 10 or higher
keep <- rowSums(counts(dds) >= 10) >= 3

vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)

rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)

sampleDists <- dist(t(assay(dds)))
sampleDists
