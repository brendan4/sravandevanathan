library(ballgown)
library(matrixStats)

gene.info <- t_data[c(9,10)]# where gene.info is subseted from 

setwd("C:\\Users\\brendan\\Documents\\sravandevanathan\\ballgown") #folder with all the files
samples.list <- list.files(path =".") # generate sample list: *** all samples included  
samples.list <- colnames(expressed.genes) # filtered out troublesome cols 


bg = ballgown(samples = samples.list, meas='FPKM') # generation of a ballgown object 
sampleNames(bg) # sanity check sample names

over200 <- exprfilter(gown = bg, cutoff= 200, meas = "FPKM") #filtering example

bg_filt = subset(bg,"rowVars(texpr(bg)) >
1",genomesubset=TRUE)

fpkm = texpr(bg_filt,meas="FPKM")
fpkm = log2(fpkm+1)
boxplot(fpkm, las = 2, ylab='log2(FPKM+1)')


