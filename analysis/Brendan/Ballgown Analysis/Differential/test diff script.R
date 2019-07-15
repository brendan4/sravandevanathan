data("full.pheno.table")
setwd("~/sravandevanathan/ballgown") #folder with all the files

remove<- full.pheno.table[which(full.pheno.table$pheno %in% c("SS", "S")),1]

#OPTION 2: filtered out troublesome cols 
samples.list <- colnames(expressed.genes)
filt.samples <- samples.list[which(!samples.list %in% remove)]

bg = ballgown(samples = filt.samples, meas='FPKM') # generation of a ballgown object 

bg_filt = subset(bg,"rowVars(texpr(bg)) >
1",genomesubset=TRUE)

full.pheno.table <- full.pheno.table[which(full.pheno.table$`colnames(expressed.genes)` %in% sampleNames(bg_filt)),]

pData(bg_filt) = data.frame(id = sampleNames(bg_filt), 
                          pheno = full.pheno.table$pheno, 
                          replicates = full.pheno.table$Replicates)
pData(bg_filt)
gene.results = stattest(bg_filt, feature="gene", covariate="pheno", meas="FPKM", getFC = TRUE)
gene.results.filt <- diff.genes.cleanup(gene.results, bg_filt, subset = TRUE)
