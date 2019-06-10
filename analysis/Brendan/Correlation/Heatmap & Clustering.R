

# hierarchical clustering of samples and heatmap of sample similarities
library(pheatmap)
pheatmap(cor(na.omit(expressed.genes), method = "spearman")) 

# Principle coordiante analysis
library(vegan)
library(rgl)
library(ape)

#simple clustering based on FPKMs
d <- cor(na.omit(expressed.genes), method="spearman")
hc <- hclust(dist(1-d))
plot.phylo(as.phylo(hc), type="p", edge.col=4, edge.width=3, show.node.label=TRUE, no.margin=TRUE)
dev.off()

# assembling table of conditions to lable PCoA plot:
# (in the chunk below, replace factor1 and factor2 with your actual factor names from myConditions table)
factor1=as.character(colData(dds)$factor1)
factor2=as.character(colData(dds)$factor2)
oneByTwo=paste(factor1,factor2,sep=".")
conditions=data.frame(cbind(factor1,factor2,oneByTwo))

# actual PCoA analysis
dds.pcoa=pcoa(vegdist(t(vsd),method="manhattan")/1000)
scores=dds.pcoa$vectors

# plotting
plot(scores[,1], scores[,2],col=as.numeric(as.factor(factor1)))
ordispider(scores,factor2,label=T)
ordiellipse(scores,factor2)

# interactive 3d plot - can rotate it by dragging mouse
radiusScale=2 # change it if the spheres come out too small or too big in the next one
plot3d(scores[,1], scores[,2],scores[,3],col=as.numeric(as.factor(factor1)),type="s",radius=radiusScale*as.numeric(as.factor(factor2)))

# formal permutation-based analysis of variance 
adonis(t(vsd)~factor1*factor2,data=conditions,method="manhattan")  




