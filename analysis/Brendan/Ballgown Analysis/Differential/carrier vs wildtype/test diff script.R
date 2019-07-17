data("full.pheno.table")
setwd("~/sravandevanathan/ballgown") #folder with all the files

#removal of diseased indivduals 
remove<- full.pheno.table[which(full.pheno.table$pheno %in% c("SS", "S")),1]

# filtered out troublesome cols from ballgown processing
samples.list <- colnames(expressed.genes)
filt.samples <- samples.list[which(!samples.list %in% remove)]

bg = ballgown(samples = filt.samples, meas='FPKM') # generation of a ballgown object 

library(matrixStats)
bg_filt = subset(bg,"rowVars(texpr(bg)) >
1",genomesubset=TRUE)

full.pheno.table <- full.pheno.table[which(full.pheno.table$`colnames(expressed.genes)` %in% sampleNames(bg_filt)),]

pData(bg_filt) = data.frame(id = sampleNames(bg_filt), 
                          pheno = full.pheno.table$pheno, 
                          replicates = full.pheno.table$Replicates)
pData(bg_filt)
gene.results = stattest(bg_filt, feature="gene", covariate="pheno", meas="FPKM", getFC = TRUE)

# optional filterings
gene.results.filt <- diff.genes.cleanup(gene.results, bg_filt, subset = TRUE)



#manually subseting data
#gene names matching
indices <- match(gene.results$id, texpr(bg_filt, 'all')$gene_id)
gene_names_for_result <- texpr(bg_filt, 'all')$gene_name[indices]
diff.results <- data.frame(geneNames = gene_names_for_result, gene.results)

diff.results = arrange(diff.results, pval)
filt.names <- subset(diff.results, diff.results$qval<0.07)

#filters outs periods
filt.set <- filter.out.genes(filt.names, gene.list = c("\\."), by.rownames = FALSE, col = 1)
diff.genes <- filter.genes(expressed.genes, filt.set$geneNames, lazy = FALSE)

#ploting
sig=which(diff.results$pval<0.05)
diff.results[,"de"] = log2(diff.results[,"fc"])
hist(diff.results[sig,"de"], breaks=50, col="seagreen", 
     xlab="log2(Fold change) Sen_DS vs Sen_WW", 
     main="Distribution of differential expression values")

abline(v=-2, col="black", lwd=2, lty=2)
abline(v=2, col="black", lwd=2, lty=2)
#optionals
legend("topleft", "Fold-change > 4", lwd=2, lty=2)


carrier <- full.pheno.table$`colnames(expressed.genes)`[which(full.pheno.table$pheno == "C")]
wildtype <- full.pheno.table$`colnames(expressed.genes)`[which(full.pheno.table$pheno == "W")]

#base r plot for noobs 
expressed.genes[,"wildtype"] <- apply(expressed.genes[,which(colnames(expressed.genes) %in% wildtype)], 1, mean)
expressed.genes[,"carrier"] <- apply(expressed.genes[,which(colnames(expressed.genes) %in% carrier)], 1, mean)
x=log2(expressed.genes[,"wildtype"]+.1)
y=log2(expressed.genes[,"carrier"]+.1)
plot(x=x, y=y, pch=16, cex=0.25, xlab="Wildtype FPKM (log2)", ylab="Carrier FPKM (log2)", main="Wildtype vs Carrier FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

#ggplot for legends
expressed.genes[,"carrier"] <- log2(expressed.genes[,"carrier"]+.1)
expressed.genes[,"wildtype"] <- log2(expressed.genes[,"wildtype"]+.1)

sig <- diff.results[which(diff.results$pval < 0.05), c(1,2)]
expressed.genes <- pretty.gene.name(expressed.genes)
sig.data <-  expressed.genes[which(expressed.genes$pretty %in% sig$geneNames), ]

#qsig 
sigq <- diff.results[which(diff.results$qval < 0.07), c(1,2)]
sigq.data <- filter.genes(expressed.genes, sigq$geneNames ,lazy = TRUE) #Maybe use lazy = FALSE

library(ggplot2)
ggplot(expressed.genes, aes(x = wildtype, y = carrier))+ 
  ggtitle("Wildtype vs Carrier FPKMs")+
  xlab("Wildtype FPKM (log2)")+
  ylab("Carrier FPKM (log2)")+
  geom_point(alpha = .2) + 
  geom_point(data = sig.data, aes(x = wildtype, y = carrier), colour = "orange", alpha = .3)+
  geom_label(data = sigq.data, 
             aes(x = wildtype, y = carrier), 
             label = sigq.data$pretty, 
             colour = "deepskyblue4", alpha =1,
             size = 1.75)
