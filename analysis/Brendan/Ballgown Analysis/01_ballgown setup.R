library(ballgown)
library(matrixStats)

setwd("~/sravandevanathan/ballgown") #folder with all the files

#OPTION 1: generate sample list: *** all samples included  
samples.list <- list.files(path =".")

#OPTION 2: filtered out troublesome cols <- RECOMENDED
samples.list <- colnames(expressed.genes)


bg = ballgown(samples = samples.list, meas='FPKM') # generation of a ballgown object 
sampleNames(bg) # sanity check sample names

#OPTIONS 1: filter out by expression 
over200 <- exprfilter(gown = bg, cutoff= 200, meas = "FPKM") #filtering example

#OPTION 2: filter out low expressed <-uses this for most general cases 
bg_filt = subset(bg,"rowVars(texpr(bg)) >
1",genomesubset=TRUE)

#filter boxplot
fpkm = texpr(bg_filt,meas="FPKM")
fpkm = log2(fpkm+1)
boxplot(fpkm, las = 2, ylab='log2(FPKM+1)')

#transcripts lengths 
full_table <- texpr(bg , 'all')
hist(full_table$length, breaks=1000, xlab="Transcript length (bp)", main="Distribution of transcript lengths", col="steelblue",  xlim = range(0,30000))
