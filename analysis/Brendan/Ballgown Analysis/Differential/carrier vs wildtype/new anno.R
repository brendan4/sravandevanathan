setwd("~/sravandevanathan/ballgown_new annotation") #folder with all the files
#please generate the pheno data which combines july and old with july pheno generation in pheno folder
full.pheno.table <- pheno.more.samples
full.pheno.table$`colnames(expressed.genes)` <- levels(droplevels(full.pheno.table$`colnames(expressed.genes)`))
full.pheno.table$`colnames(expressed.genes)`[33:36] <- paste("121317_", full.pheno.table$`colnames(expressed.genes)`[33:36], sep = "")

#removal of diseased indivduals 
remove<- full.pheno.table[which(full.pheno.table$pheno %in% c("SS", "S")),1]

# filtered out troublesome cols from ballgown processing
samples.list <- colnames(expressed.genes.GEN)
filt.samples <- samples.list[which(!samples.list %in% remove)]
filt.samples[29:32] <- paste("121317_", filt.samples[29:32], sep = "")

bg = ballgown(samples = filt.samples, meas='all') # generation of a ballgown object 

#filtering of data
bg_filt <- exprfilter(gown = bg, cutoff= .75, meas = "FPKM") 

#pheno data prep
full.pheno.table <- full.pheno.table[which(full.pheno.table$`colnames(expressed.genes)` %in% sampleNames(bg_filt)),]
pData(bg_filt) = data.frame(id = sampleNames(bg_filt), 
                            pheno = full.pheno.table$pheno, 
                            replicates = full.pheno.table$Replicates)
pData(bg_filt)

#diff expression
gene.results = stattest(bg_filt, feature="transcript", covariate="pheno", meas="FPKM", getFC = TRUE)

diff.results = data.frame(geneNames=ballgown::geneNames(bg_filt), geneIDs=ballgown::geneIDs(bg_filt), gene.results)

# optional result filterings
gene.results.filt <- diff.genes.cleanup(gene.results, bg_filt, subset = TRUE)

#manually subseting data 
#gene names matching

indices <- match(gene.results$id, texpr(bg_filt, 'all')$gene_id)
gene_names_for_result <- texpr(bg_filt, 'all')$gene_name[indices]
diff.results <- data.frame(geneNames = gene_names_for_result, gene.results)

#signficance filtering 
diff.results = arrange(diff.results, pval)
filt.names <- subset(diff.results, diff.results$qval<0.07)
filt.set <- filter.out.genes(filt.names, gene.list = c("\\."), by.rownames = FALSE, col = 1)
diff.genes <- filter.genes(expressed.genes.GEN, filt.set$geneNames, lazy = FALSE)

#ploting hist dist fc 
sig <- which(diff.results$pval<0.05)
diff.results[,"de"] = log2(diff.results[,"fc"])
hist(diff.results[sig,"de"], breaks=50, col="seagreen", 
     xlab="log2(Fold change) Sen_DS vs Sen_WW", 
     main="Distribution of differential expression values")

abline(v=-2, col="black", lwd=2, lty=2)
abline(v=2, col="black", lwd=2, lty=2)
#optionals
legend("topleft", "Fold-change > 4", lwd=2, lty=2)

#plot for mean wildtype vs carrier: data prep
carrier <- full.pheno.table$`colnames(expressed.genes)`[which(full.pheno.table$pheno == "C")]
wildtype <- full.pheno.table$`colnames(expressed.genes)`[which(full.pheno.table$pheno == "W")]
expressed.genes.GEN[,"wildtype"] <- apply(expressed.genes.GEN[,which(colnames(expressed.genes.GEN) %in% wildtype)], 1, mean)
expressed.genes.GEN[,"carrier"] <- apply(expressed.genes.GEN[,which(colnames(expressed.genes.GEN) %in% carrier)], 1, mean)

#base r plot for noobs 
x=log2(expressed.genes.GEN[,"wildtype"]+.1)
y=log2(expressed.genes.GEN[,"carrier"]+.1)
plot(x=x, y=y, pch=16, cex=0.25, xlab="Wildtype FPKM (log2)", ylab="Carrier FPKM (log2)", main="Wildtype vs Carrier FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

#data prep for ggplots 
sig <- diff.results[which(diff.results$pval < 0.05), c(1,2)]
expressed.genes.GEN <- pretty.gene.name(expressed.genes.GEN)
sig.data <-  expressed.genes.GEN[which(expressed.genes.GEN$pretty %in% sig$geneNames), ]

#qsig 
sigq <- diff.results[which(diff.results$qval < 0.07), c(1,2)]
sigq.data <- filter.genes(expressed.genes.GEN, sigq$geneNames ,lazy = TRUE) #Maybe use lazy = FALSE

#ggplot for legends
expressed.genes.GEN[,"carrier"] <- log2(expressed.genes.GEN[,"carrier"]+.1)
expressed.genes.GEN[,"wildtype"] <- log2(expressed.genes.GEN[,"wildtype"]+.1)


library(ggplot2)
ggplot(expressed.genes.GEN, aes(x = wildtype, y = carrier))+ 
  ggtitle("Wildtype vs Carrier FPKMs")+
  xlab("Wildtype FPKM (log2)")+
  ylab("Carrier FPKM (log2)")+
  geom_point(alpha = .2) + 
  geom_point(data = sig.data, aes(x = log2(wildtype), y = log2(carrier)), colour = "orange", alpha = .3)+
  geom_label(data = sigq.data, 
             aes(x = log2(wildtype), y = log2(carrier)), 
             label = sigq.data$pretty, 
             colour = "deepskyblue4", alpha =1,
             size = 2)

#mean expressions for all genes 
expressed.genes.GEN <- expressed.genes.GEN[,-which(colnames(expressed.genes.GEN) %in% c("means", "pretty"))]
expressed.genes.GEN[,"means"] <- apply(expressed.genes.GEN, 1, mean)

# duplicate handling and matching 
exp.filter <- pretty.gene.name(expressed.genes.GEN)
exp.filter <- exp.filter[!duplicated(exp.filter$pretty),]
diff.filter <- diff.results[!duplicated(diff.results$geneNames),]
exp.filter <- exp.filter[which(exp.filter$pretty %in% diff.results$geneNames),]
diff.filter <- diff.results[which(diff.results$geneNames %in% exp.filter$pretty),]


#sig data prep
sig <- diff.filter[which(diff.filter$pval < 0.001), c(1,2,4)]
expressed.genes.GEN <- pretty.gene.name(expressed.genes.GEN)
sig.data <-  exp.filter[which(exp.filter$pretty %in% sig$geneNames), ]
sig.data <- sig.data[order(sig.data$pretty),]
sig <- sig[order(sig$geneNames), ]
sig.data$fc <- sig$fc

#qsig 
sigq <- diff.filter[which(diff.filter$qval < 0.07), c(1,2,4)]
sigq.data <- filter.genes(exp.filter, sigq$geneNames ,lazy = TRUE) #Maybe use lazy = FALSE
sigq.data <- sigq.data[order(sigq.data$pretty),]
sigq <- sigq[order(sigq$geneNames), ]
sigq.data$fc <- sigq$fc
sigq.data <- sigq.data[-grep("RPL11P3", rownames(sigq.data)),]

diff.filter <- diff.filter[!(duplicated(diff.filter$geneNames)),]

ggplot(exp.filter, aes(x = log2(exp.filter$means), y = log2(diff.filter$fc)))+
  xlab("Average log2FPKM")+
  ylab("log2FC")+
  geom_point()+ 
  geom_point(data = sig.data, aes(x = log2(means), y = log2(fc)), 
             colour = "orange", 
             alpha = .8)+
  geom_label(data = sigq.data, 
             aes(x = log2(means), y = log2(fc)), 
             label = sigq.data$pretty, 
             colour = "deepskyblue4", alpha =1,
             size = 2)

# scale not log
ggplot(exp.filter, aes(x = exp.filter$means, y = diff.filter$fc))+
  xlab("Average FPKM")+
  ylab("FC")+
  scale_x_continuous(trans = 'log2') +
  scale_y_continuous(trans = 'log2') +
  geom_point()+ 
  geom_point(data = sig.data, aes(x = means, y = fc), 
             colour = "skyblue3", 
             alpha = .8)+
  geom_label(data = sigq.data, 
             aes(x = means, y = fc), 
             label = sigq.data$pretty, 
             colour = "navy", alpha =1,
             size = 2)

#valcano plots
#data subsets 
sigp <- diff.results[which(diff.results$pval < 0.05),]
sigq <- diff.results[which(diff.results$qval < 0.07),]

ggplot(data = diff.results, aes(x= log2(fc), y= -log2(pval)))+
  xlab("log2FC")+ 
  ylab("-log2pval")+
  geom_point()+
  geom_point(data = sigp, aes(x= log2(fc), y= -log2(pval)), 
             color = "skyblue3")+
  geom_label(data = sigq, aes(x= log2(fc), y = -log2(pval)), 
             label = sigq$geneNames, 
             color = "navy", size = 2.5, position= "jitter")

library(Glimma)
diff.results[grep("BAG1", diff.results$geneNames),]
expressed.genes.GEN[grep("BAG1", rownames(expressed.genes.GEN)),]


results_transcripts = data.frame(geneNames=ballgown::geneNames(bg_filt), geneIDs=ballgown::geneIDs(bg_filt), gene.results)
#valcano plots
#data subsets 
sigp <- results_transcripts[which(results_transcripts$pval < 0.05),]
sigq <- results_transcripts[which(results_transcripts$qval < 0.07),]

ggplot(data = results_transcripts, aes(x= log2(fc), y= -log2(pval)))+
  xlab("log2FC")+ 
  ylab("-log2pval")+
  geom_point()+
  geom_point(data = sigp, aes(x= log2(fc), y= -log2(pval)), 
             color = "skyblue3")+
  geom_label(data = sigq, aes(x= log2(fc), y = -log2(pval)), 
             label = sigq$geneNames, 
             color = "navy", size = 2.5, position= "jitter")