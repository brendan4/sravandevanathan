
min_nonzero = 1
gene.list <- c("MDM2","RPL11", "TP53", "GATA1", "PML", "MYC", "CDKN2A")
# two options: expressed genes and expressed trans datasets
filtered.data <- filter.genes(expressed.genes,gene.list) # expressed genes dataset 
filtered.data <- filter.genes(expressed.trans,gene.list) # expressed trans dataset 


#24 interesting protiens from gene.list
par(mar=c(7,4,4,2)+0.1) 
png(filename='Distribution of FPKMs.png', width=800, height=750)
boxplot(log2(filtered.data+min_nonzero), 
        names=colnames(filtered.data), las=2, ylab="log2(FPKM)", 
        main="Distribution of FPKMs for all libraries")
graphics.off()

#distribution in Ribo protiens 
par(mar=c(7,4,4,2)+0.1) 
png(filename='Distribution of Ribo  FPKMs.png', width=800, height=750)
boxplot(log2(Ribo.filter+min_nonzero), 
        names=colnames(filtered.data), las=2, ylab="log2(FPKM)", 
        main="Distribution of Ribo FPKMs")
graphics.off()

sub <- rownames(na.omit(filtered.data[which(filtered.data > 100),])) # RPL11 highest abundance in trans
RPL11 <- expressed.trans[grep("^RPL11",rownames(expressed.trans)),] # all RPL11 in trans 

#text(Ribo.filter, Ribo.filter,RPL11, cex=0.6, pos=4, col="red")

#a look at myc 
filt <- expressed.genes[grep("^MYC", rownames(expressed.genes)),]
filt$var <- apply(log2(filt +0.01),1, var)
filt <- pretty.gene.name(filt)
rownames(filt) <- filt$pretty
filt.p <- filt[,-which(colnames(filt) %in% c("pretty"))]
ggplot(filt.p, aes(x = rownames(filt.p) , y = var)) + geom_bar(stat="identity", fill="tomato3") 

# a look at RPL11
filt <- expressed.genes[grep("^RPL11", rownames(expressed.genes)),]
filt$var <- apply(log2(filt +0.01),1, var)
filt <- pretty.gene.name(filt)
rownames(filt) <- filt$pretty
filt.p <- filt[,-which(colnames(filt) %in% c("pretty"))]

#PML 
filt <- expressed.genes[grep("^PML", rownames(expressed.genes)),]
filt$var <- apply(log2(filt +0.01),1, var)
filt <- pretty.gene.name(filt)
rownames(filt) <- filt$pretty
filt.p <- filt[,-which(colnames(filt) %in% c("pretty"))]


#a look at all genes in gene.list 

filtered.data$var <- apply(log10(filtered.data),1, var)
filt <- pretty.gene.name(filtered.data)
rownames(filt) <- filt$pretty
filt.p <- filt[,-which(colnames(filt) %in% c("pretty"))]
filt.p <- na.omit(filt.p)
ggplot(filt.p, aes(x = rownames(filt.p) , y = var)) + geom_bar(stat="identity", fill="tomato3") 
colSum <- as.data.frame(colSums(filt.p))
colSum <- t(colSum)
filt.p <- rbind(filt.p, colSum)


# a look a hemo protiens 
gene.list <- c('HBA', "HBB", "HBG1", "HBG2", "HBE", "HBD", "SLC4A1", "SNCA", "BPGM")
filtered.data <- filter.genes(expressed.genes, gene.list = gene.list)
filtered.data$var <- apply(log2(filtered.data +0.1),1, var)
filt <- pretty.gene.name(filtered.data)
rownames(filt) <- filt$pretty
filt.p <- filt[,-which(colnames(filt) %in% c("pretty"))]
ggplot(filt.p, aes(x = rownames(filt.p) , y = var)) + geom_bar(stat="identity", fill="tomato3") 

# HB
filtered.data <- filter.genes(expressed.genes, "HB")
filtered.data$var <- apply(log2(filtered.data +0.1),1, var)
filt <- pretty.gene.name(filtered.data)
rownames(filt) <- filt$pretty
filt.p <- filt[,-which(colnames(filt) %in% c("pretty"))]
ggplot(filt.p, aes(x = rownames(filt.p) , y = var)) + geom_bar(stat="identity", fill="tomato3") 

# H2BFXP
filtered.data <- filter.genes(expressed.genes, "H2BFXP")
filtered.data$var <- apply(log2(filtered.data +0.1),1, var)
filt <- pretty.gene.name(filtered.data)
rownames(filt) <- filt$pretty
filt.p <- filt[,-which(colnames(filt) %in% c("pretty"))]


#FNDC4
filtered.data <- filter.genes(expressed.genes, "FNDC4")
filtered.data$var <- apply(log2(filtered.data +0.1),1, var)
filt <- pretty.gene.name(filtered.data)
rownames(filt) <- filt$pretty
filt.p <- filt[,-which(colnames(filt) %in% c("pretty"))]



werid <- c("NPHP1", "MORN4","POTEI", "GPRC5B")
#immune proteins 
filtered.data <- filter.genes(expressed.genes, c("JCHAIN","MX1","IGLL5","RSAD2","CMPK2", "IFI44L", "MZB1", "IFIT1", "OAS3"))
filtered.data$var <- apply(log2(filtered.data +0.1),1, var)
filt <- pretty.gene.name(filtered.data)
rownames(filt) <- filt$pretty
filt.p <- filt[,-which(colnames(filt) %in% c("pretty"))]
ggplot(filt.p, aes(x = rownames(filt.p) , y = var)) + geom_bar(stat="identity", fill="tomato3") 

# AHSP, FAM132B, HEMGN, and TRIM10: from the GATA1 diff DBA paper : blood cel formation 
filtered.data <- filter.genes(expressed.genes, c("AHSP", "FAM132B", "HEMGN", "TRIM10"))
filtered.data$var <- apply(log2(filtered.data +0.1),1, var)
filt <- pretty.gene.name(filtered.data)
rownames(filt) <- filt$pretty
filt.p <- filt[,-which(colnames(filt) %in% c("pretty"))]
ggplot(filt.p, aes(x = rownames(filt.p) , y = var)) + geom_bar(stat="identity", fill="tomato3") 

# heme proteins ALAS2, FECH, CPOX, PPOX, and UROS.: GATA1 paper

filtered.data <- filter.genes(expressed.genes, c("ALAS2", "FECH", "CPOX", "PPOX", "UROS"))
filtered.data$var <- apply(log2(filtered.data +0.1),1, var)
filt <- pretty.gene.name(filtered.data)
rownames(filt) <- filt$pretty
filt.p <- filt[,-which(colnames(filt) %in% c("pretty"))]
ggplot(filt.p, aes(x = rownames(filt.p) , y = var)) + geom_bar(stat="identity", fill="tomato3") 

#IL8, IL1R1, CXCR4, ICAM3, MPO, TNFSF10, and TLR4 genes with IL6, TNF,

filtered.data <- filter.genes(expressed.genes, c("IL8", "IL1R1", "CXCR4", "ICAM3", "MPO", "TLR4", "IL6", "TNF"))
filtered.data$var <- apply(log2(filtered.data +0.1),1, var)
filt <- pretty.gene.name(filtered.data)
rownames(filt) <- filt$pretty
filt.p <- filt[,-which(colnames(filt) %in% c("pretty"))]
ggplot(filt.p, aes(x = rownames(filt.p) , y = var)) + geom_bar(stat="identity", fill="tomato3") 



#RNR1 
filtered.data <- filter.genes(expressed.genes, c("RNA"))
filtered.data$var <- apply(log2(filtered.data +0.1),1, var)
filt <- pretty.gene.name(filtered.data)
rownames(filt) <- filt$pretty
filt.p <- filt[,-which(colnames(filt) %in% c("pretty"))]
ggplot(filt.p, aes(x = rownames(filt.p) , y = var)) + geom_bar(stat="identity", fill="tomato3") 


filtered.data <- Var.samples(expressed.genes, gene.list = c("RPL11"), pretty.names = TRUE, graph = TRUE)
