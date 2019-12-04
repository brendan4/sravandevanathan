library("GenomicAlignments")
library("BiocParallel")
library("GenomicFeatures")
library("Rsamtools")

# set working directory: indir should contain bam files
indir <- "F:/dba/DBA_121317" #CHANGE TO BAM files destination
ids <- list.files(setwd("~/sravandevanathan/ballgown"))

# idendifying bam files needed
filenames <- file.path(indir, paste0(pheno$id, ".bam"))
file.exists(filenames)

# interface for BAM files: yieldSize = process 2 million reads at a time
bamfiles <- BamFileList(filenames, yieldSize=2000000)
seqinfo(bamfiles[1])

# defining genome model
gtffile <- file.path(indir,"Gencode.gtf") # Assumes GTF is in the same file as bam files and has the name Gencode.gtf
txdb <- makeTxDbFromGFF(gtffile, format = "gtf", circ_seqs = character())

# GRangesList of all the exons grouped by gene
ebg <- exonsBy(txdb, by="gene")

# creating summarized experiments objects with counts
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )

hold <- assay(se)
dim(se)
