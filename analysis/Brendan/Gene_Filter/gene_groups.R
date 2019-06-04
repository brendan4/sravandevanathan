
#gene groups 
gata1 <- expressed.genes[pmatch("GATA1",expressed.genes$Gene.Names),]
p53 <- expressed.genes[grep("^TP53", expressed.genes$Gene.Names),]
TMA7 <- expressed.genes[pmatch("TMA7",expressed.genes$Gene.Names),]
apaf1 <- expressed.genes[pmatch("APAF1",expressed.genes$Gene.Names),]
#RPS19 <- expressed.genes[grep("^RPS19",expressed.genes$Gene.Names),]
#ADA <- expressed.genes[grep("^ADA",expressed.genes$Gene.Names),]
RPL11 <- expressed.genes[grep("^RPL11",rownames(expressed.genes)),]
MDM2 <- expressed.genes[grep("^MDM2",rownames(expressed.genes), ignore.case = TRUE, value = TRUE), ]

min_nonzero = 1
gene.list <- c("MDM2","RPL11", "APAf1", "TP53", "GATA1", "PML")
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
       