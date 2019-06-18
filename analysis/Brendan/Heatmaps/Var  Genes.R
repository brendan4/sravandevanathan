library(gplots)
library(ComplexHeatmap)


#estimate the variance for each row in the logcounts matrix
var_genes <- apply(log2(na.omit(expressed.genes)+0.1), 1, var)
head(var_genes)

# Get the gene names for the top 100 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:100]
head(select_var)

# Subset logcounts matrix
highly_variable <- expressed.genes[select_var,]
dim(highly_variable)

#OPTIONAL look
highly_variable$var<- apply(log2(highly_variable+0.1), 1, var)
head(var_genes)

###################SKIP: if don't want sex related genes filtered out######################
#sex insight
sex <- t(sex)
sex <- highly_variable[grep("^XIST", rownames(highly_variable)),]
sex[,which(sex > 10)] = "F"
sex[,which(sex < 10)] = "M"
rownames(sex) <- "Gender"
sex <- t(sex)

#bg_filt diff expression with sex as covariate: uses fitler ot func in gene filter/ tools
setwd("C:/Users/brendan/Documents/sravandevanathan/analysis/Brendan/Working Dataset/02_Filtered Data")
filt.names <- read.delim("sex_related.tab", sep= "")
sex <- read.delim("sample_gender.tab", sep= "")

sex.filt <- filter.out.genes(na.omit(expressed.genes), filt.names$geneNames)
#OPTIONAL: filters out LOC and unnamed genes
sex.filt <- filter.out.genes(sex.filt, c("LOC","-"))

# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(log2(sex.filt), 1, var)
head(var_genes)

# Get the gene names for the top 100 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:100]
head(select_var)

# Subset logcounts matrix
highly_variable <- expressed.genes[select_var,]
dim(highly_variable)

#OPTIONAL: make gene names pretty: load pretty func in brendan/tools 
highly_variable <- pretty.gene.name(highly_variable) # make the pretty col 
if (duplicated(highly_variable$pretty == TRUE)) {
  highly_variable <- highly_variable[-which(duplicated(highly_variable$pretty)),] # drop the duplicated 
}
rownames(highly_variable) <- highly_variable$pretty # make gene names pretty
highly_variable <- highly_variable[,-which(colnames(highly_variable) %in% c("pretty"))] # drop pretty col



par(mar=c(7,4,4,2)+0.1) 
png(filename='100 Variable Genes - sex filtered- no name filt.png', width=800, height=750)

# use t() for k means clustering of sample no genes 
Heatmap(as.matrix(log2(highly_variable+0.001)),
        column_names_side = "bottom",
        row_names_side = "left",
        row_hclust_side = "left",
        row_names_gp=gpar(cex=0.6),
        row_hclust_width = unit(3, "cm"),
        clustering_distance_rows ="euclidean",
        clustering_method_rows = "centroid",
        km=3) # number of clusters you want
graphics.off()


#simple PCA plots: call PCA func 


