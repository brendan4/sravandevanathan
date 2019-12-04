data("full.pheno.table")
data(pheno.basic)

#ball gown setup
test <- colnames(expressed.genes)
test <- test[-which( test %in% c("L3_TTAGGC", "L3_TAGCTT", "L3_GGCTAC", "L3_GATCAG"))]
samples.list <- test

setwd("C:\\Users\\brendan\\Documents\\sravandevanathan\\ballgown") #folder with all the files
bg = ballgown(samples = samples.list, meas='FPKM') # generation of a ballgown object 

#filter out low expressed
bg_filt = subset(bg,"rowVars(texpr(bg)) >
                 1",genomesubset=TRUE)

#filter means
over1 = exprfilter(bg, cutoff=20) 

full.pheno.table <- full.pheno.table[which(full.pheno.table$`colnames(expressed.genes)` %in% sampleNames(bg_filt)),]
pheno.basic <- pheno.basic[(pheno.basic$`colnames(expressed.genes)` %in% sampleNames(bg_filt)),]

full.pheno.table[which(full.pheno.table$pheno == "S"),"pheno"] <- "SS"

pData(over1) = data.frame(id = sampleNames(over1), 
                            pheno = full.pheno.table$pheno, 
                            replicates = full.pheno.table$Replicates)
pData(bg_filt)


gene.results = stattest(over1, feature="gene", covariate="pheno", meas="FPKM", adjustvars = c("replicates"))
gene.results.filt <- diff.genes.cleanup(gene.results, bg_filt, subset = TRUE )
