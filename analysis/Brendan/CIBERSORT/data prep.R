# cybersort data table prep
expressed.genes.cyber <- na.omit(expressed.genes)
#OPTIONAL: remove _ and numbers from gene names 
expressed.genes.cyber <- pretty.gene.name(expressed.genes.cyber, as.row.names = TRUE, remove.dups = TRUE)
expressed.genes.cyber$Gene<- rownames(expressed.genes.cyber)
#rownames editing
expressed.genes.cyber <- expressed.genes.cyber[,c(35,1:34)]
rownames(expressed.genes.cyber) <- c(seq(1,nrow(expressed.genes.cyber)))
# col names editing 
expressed.genes.cyber[(nrow(expressed.genes.cyber)+1),] <- colnames(expressed.genes.cyber)
expressed.genes.cyber <- expressed.genes.cyber[c(nrow(expressed.genes.cyber),1:(nrow(expressed.genes.cyber)-1)),]
colnames(expressed.genes.cyber) <- c(seq(1,ncol(expressed.genes.cyber)))
rownames(expressed.genes.cyber) <- c(seq(1,nrow(expressed.genes.cyber)))
#OPTIONAL: resets row names
rownames(expressed.genes.cyber) <- expressed.genes.cyber[,1]
# saving 
expressed.genes.cyber <- as.data.frame(expressed.genes.cyber)
write.csv(expressed.genes.cyber,"expressed_genes_cyber.csv",row.names = FALSE, col.names = FALSE)
new <- read.table("expressed_genes_cyber.txt")

# removes - from row names 
expressed.genes.cyber[grepl("^-", expressed.genes.cyber$`1`),1] <-  substring(expressed.genes.cyber[grepl("^-", expressed.genes.cyber$`1`),1], 
                                                                              first = 2) 

#saving test
write.csv(ExampleMixtures.GEPs,"expressed_genes_cyber.csv",row.names = FALSE, col.names = NA)
new <-read.csv("expressed_genes_cyber.csv")
