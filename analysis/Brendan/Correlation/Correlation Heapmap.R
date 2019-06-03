# this script uses the Corelation Matrix.R script: load cor.table function from that script

cor.matrix.genes <- cor.table(na.omit(expressed.genes))
cor.matrix.trans <- cor.table(na.omit(expressed.trans))

library(gplots)
library(ggcorrplot)

ggcorrplot(cor.matrix.genes, hc.order = TRUE, type = "lower", colors = c("blue","white", "red"))

#image saving options
par(mar=c(7,4,4,2)+0.1) 
png(filename='.png', width=800, height=750)

graphics.off()
