URL <- "https://www.genenames.org/cgi-bin/genegroup/download?id=1054&type=branch"
Ribo.protiens <- read.table(URL, sep = '\t',header = TRUE) # importing 
Ribo.protiens <- Ribo.protiens[,c(2,7)] # sybolms and alternative names 
Ribo.names <- t(Ribo.protiens) # transpose
rm(Ribo.protiens) # removes rib.protien from the workspace
Ribo.names <- Ribo.names[1,] # selection of symbols 
Ribo.names <- Filter(function(x) !any(grepl("^MRP", x)), Ribo.names) # drop mitochondrial genes: optional prob. needed to limit gene number
Ribo.names <- Filter(function(x) any(grepl("^RPL", x)), Ribo.names)  # only large: optional 
Ribo.names <- Filter(function(x) any(grepl("^RPS", x)), Ribo.names)  # only small: optional

Ribo.filter <- data.frame() # intialize data frame

#selection of genes related to ribosomes
for (name in Ribo.names){
  print(paste("^",name, sep= "")) # prints out pattern
  match<- expressed.genes[grep(paste("^", name, sep = ""),rownames(expressed.genes), ignore.case = TRUE, value = TRUE),]
  Ribo.filter <- rbind(Ribo.filter, match) # merges to preivous 
}

library(ggplot2)
library(devtools)


#princple component for expressed.genes
genes.PCA <- prcomp(scale(t(na.omit(Ribo.filter)))) # t for gene comp
summary(genes.PCA)

#basic plot of PCA1 and PCA 2
plot(genes.PCA$x[,1], genes.PCA$x[,2],
     main = "PCA: Genes - log scale",
     xlab = paste("PC1 - ", genes.PCA.var.per[1], "%", sep = ""), 
     ylab = paste("PC2 - ", genes.PCA.var.per[2], "%", sep = ""))

#precent variance 
genes.PCA.var <- genes.PCA$sdev^2
genes.PCA.var.per <- round(genes.PCA.var/sum(genes.PCA.var)*100, 1)

#scree plots of variation 
barplot(genes.PCA.var.per, main = "Scree Plot: Genes", 
        xlab = paste("PC1 - ", genes.PCA.var.per[1], "%", sep = ""), 
        ylab = paste("PC2 - ", genes.PCA.var.per[2], "%", sep = ""))
       

# reformate of PCA data into a data frame for ggplot usage
genes.PCA.data <- data.frame(Sample = rownames(genes.PCA$x),
                             x = genes.PCA$x[,1],
                             y = genes.PCA$x[,2])

#ggplot of PCA data

ggplot(data = genes.PCA.data, aes(x = x, y = y, label = Sample))+
  geom_text(size = 2)+
  xlab(paste("PC1 - ", genes.PCA.var.per[1], "%", sep = ""))+
  ylab(paste("PC2 - ", genes.PCA.var.per[2], "%", sep = ""))+
  ggtitle("PCA: Expressed Genes")

#gene_PCA_plot <- ggbiplot(genes.PCA)

# 100 genes the influence the PCA the greatest (either pos of neg)
gene_score_ranked <- sort(abs(genes.PCA$rotation[,1]))
gene_top_hun <- names(gene_score_ranked[1:100]) # reduce for this data type
genes.PCA$rotation[gene_top_hun, 1] # push to left on x axis (-) or right on x axis (+)

