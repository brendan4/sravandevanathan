bg_filt = subset(bg,"rowVars(texpr(bg)) >
                 1",genomesubset=TRUE)

fpkm = texpr(bg_filt,meas="FPKM")
fpkm = log2(fpkm+1)
boxplot(fpkm, las = 2, ylab='log2(FPKM+1)')


