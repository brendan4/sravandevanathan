# Classical MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name

library(ggplot2)

#distance - Euclidean 
#d_genes <- dist(scale(t(expressed.genes),center = TRUE, scale = TRUE), method = "euclidean") # for scaled values
d_genes <- dist(t(expressed.genes), method = "euclidean")

#MDS with euclidean distance
mds_genes <- cmdscale(d_genes, eig = TRUE, x.ret = TRUE)

#variation percentage 
mds.var.per_genes <- round(mds_genes$eig/sum(mds_genes$eig)*100,1)

#organization of data for ggplot use
mds.values_genes  <- mds_genes$points
mds.data_genes <- data.frame(Sample = rownames(mds.values_genes),
                             X = mds.values_genes[,1],
                             Y = mds.values_genes[,2])

#plot of MDS
ggplot(data = mds.data_genes, aes(x=X, y=Y, label = Sample))+
  geom_text(size = 2.5)+
  theme_bw()+
  xlab(paste("MDS1 - ", mds.var.per_genes[1], "%", sep = ""))+
  ylab(paste("MD2 - ", mds.var.per_genes[2], "%", sep = ""))+
  ggtitle("MDS using Euclidean distance - Genes")


#seconds mds plot avg abs value of the log fold change

log2.data.matrix.genes <- log2(na.omit(expressed.genes)+0.1) # 0.1 for zeros, NAs will mutate values
#creation of a distance matrix: with zeros 
log2.distance.matrix.genes <- matrix(0,
                                     nrow=ncol(log2.data.matrix.genes),
                                     ncol=ncol(log2.data.matrix.genes),
                                     dimnames = list(colnames(log2.data.matrix.genes),
                                                     colnames(log2.data.matrix.genes)))
# fill the table: compution of distance
for ( i in 1:ncol(log2.distance.matrix.genes)){
  for(j in 1:i){
    print(mean(abs(log2.data.matrix.genes[,i]-log2.data.matrix.genes[,j])))
    log2.distance.matrix.genes[i, j] <- mean(abs(log2.data.matrix.genes[,i]-log2.data.matrix.genes[,j]))
  }
}

mds.stuff_genes <- cmdscale(as.dist(log2.distance.matrix.genes),
                            eig = TRUE,
                            x.ret = TRUE)
mds.var.per_genes2 <- round(mds.stuff_genes$eig/sum(mds.stuff_genes$eig)*100,1)

mds.values_genes2  <- mds.stuff_genes$points
mds.data_genes2 <- data.frame(Sample = rownames(mds.values_genes2),
                                X = mds.values_genes2[,1],
                                Y = mds.values_genes2[,2])

ggplot(data = mds.data_genes2, aes(x=X, y=Y, label = Sample))+
  geom_text(size = 2.5)+
  theme_bw()+
  xlab(paste("MDS1 - ", mds.var.per_genes2[1], "%", sep = ""))+
  ylab(paste("MD2 - ", mds.var.per_genes2[2], "%", sep = ""))+
  ggtitle("MDS - Genes")
