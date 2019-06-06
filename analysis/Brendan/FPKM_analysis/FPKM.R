expressed.genes <- read.table("expressed.genes.tab") # will read the table 
expressed.trans<- read.table("expressed.trans.ctab") # will read the table 

#********Gene abundance file **********
#Set the minimum non-zero FPKM values for use later.
min_nonzero = 1 

#**DONT RUN MORE THAN ONCE**
colnames(expressed.genes) <- substring(colnames(expressed.genes), first = 6) # removes FPKM. from row names 

#log2 distribution in FPKM values 
par(mar=c(7,4,4,2)+0.1) 
png(filename='Distribution of FPKMs.png', width=800, height=750)
boxplot(log2(expressed.genes+min_nonzero), 
        names=colnames(expressed.genes), las=2, ylab="log2(FPKM)", 
        main="Distribution of FPKMs for all libraries")
graphics.off()

#plotting a suspected replicate pairs from PCA and MDS 
expressed.genesNA <- na.omit(expressed.genes)
x = expressed.genesNA[,"L6_ACAGTG"]
y = expressed.genesNA[,"L2_GCCAAT"]

plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), 
     pch=16, col="blue", cex=0.25, 
     xlab="FPKM (L_ACAGTG, Replicate 1)", 
     ylab="FPKM (L2_GCCAAT, Replicate 2)", 
     main="Comparison of expression values for a pair of replicates")

abline(a=0,b=1)
rs=cor(x,y)^2
legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")

#heat map of suspected replicates 
colors = colorRampPalette(c("white", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
smoothScatter(x=log2(x+min_nonzero), y=log2(y+min_nonzero), 
              xlab="FPKM (L2_CAGATC, Replicate 1)", 
              ylab="FPKM (L6_GCCAAT, Replicate 2)", 
              main="Comparison of expression values for a pair of replicates", 
              colramp=colors, nbin=200)


#corelation distance
expressed.genes[,"sum"]=apply(expressed.genes, 1, sum) # generates max sum across samples for each transcript
max(na.omit(expressed.genes$sum)) # max sum (across samples) for FPKM values 
i = which(expressed.genes[,"sum"] > 5) # filter out anything lower than 5 FPKM expression across the board
r=cor(expressed.genes[i,], use="pairwise.complete.obs", method="pearson") # corelation between all pairs of data
expressed.genes <- expressed.genes[,-which(colnames(expressed.genes) %in% c("sum"))] # drops sum col

# generations of an MDS plot 
d=1-r
mds=cmdscale(d, k=2, eig=TRUE)
par(mfrow=c(1,1))
plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (all non-zero genes) for all libraries", xlim=c(-0.15,0.15), ylim=c(-0.15,0.15))
points(mds$points[,1], mds$points[,2], col="grey", cex=2, pch=16)
text(mds$points[,1], mds$points[,2], colnames(expressed.genes))


#***********Transcripts file***********
#Set the minimum non-zero FPKM values for use later.
min_nonzero = 1 

#**WARNING**DONT RUN MORE THAN ONCE**
colnames(expressed.trans) <- substring(colnames(expressed.trans), first = 6) # removes FPKM. from row names 

#log2 distribution in FPKM values 
par(mar=c(7,4,4,2)+0.1) 
png(filename='Distribution of FPKMs.png', width=800, height=750)
boxplot(log2(expressed.genes+min_nonzero), 
        names=colnames(expressed.genes), las=2, ylab="log2(FPKM)", 
        main="Distribution of FPKMs for all libraries")
graphics.off()

#plotting a suspected replicate pairs from PCA and MDS 
expressed.transNA <- na.omit(expressed.trans)
x = expressed.trans[,"L2_CAGATC"]
y = expressed.trans[,"L6_GCCAAT"]

plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), 
     pch=16, col="blue", cex=0.25, 
     xlab="FPKM (L2_CAGATC, Replicate 1)", 
     ylab="FPKM (L6_GCCAAT, Replicate 2)", 
     main="Comparison of expression values for a pair of replicates")

abline(a=0,b=1)
rs=cor(x,y)^2
legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")

#heat map of suspected replicates 
colors = colorRampPalette(c("white", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
smoothScatter(x=log2(x+min_nonzero), y=log2(y+min_nonzero), 
              xlab="FPKM (L2_CAGATC, Replicate 1)", 
              ylab="FPKM (L6_GCCAAT, Replicate 2)", 
              main="Comparison of expression values for a pair of replicates", 
              colramp=colors, nbin=200)


#corelation distance
expressed.trans[,"sum"]=apply(expressed.trans, 1, sum) # generates max sum across samples for each transcript
max(na.omit(expressed.trans$sum)) # max sum (across samples) for FPKM values 
i = which(expressed.trans[,"sum"] > 5) # filter out anything lower than 5 FPKM expression across the board
r=cor(expressed.trans[i,], use="pairwise.complete.obs", method="pearson") # corelation between all pairs of data
expressed.trans <- expressed.trans[,-which(colnames(expressed.trans) %in% c("sum"))] # drops sum col

# generations of an MDS plot 
d=1-r
mds=cmdscale(d, k=2, eig=TRUE)
par(mfrow=c(1,1))
plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (all non-zero genes) for all libraries", xlim=c(-0.15,0.15), ylim=c(-0.15,0.15))
points(mds$points[,1], mds$points[,2], col="grey", cex=2, pch=16)
text(mds$points[,1], mds$points[,2], colnames(expressed.trans))

rm(d,r,mds,i,rs,x,y,colors)
