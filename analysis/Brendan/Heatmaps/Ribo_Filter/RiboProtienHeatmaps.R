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

#selection of genes realted to ribosomes
for (name in Ribo.names){
  print(paste("^",name, sep= "")) # prints out pattern
  match<- expressed.genes[grep(paste("^", name, sep = ""),rownames(expressed.genes), ignore.case = TRUE, value = TRUE),]
  Ribo.filter <- rbind(Ribo.filter, match) # merges to preivous 
}

# color and packages 
require("RColorBrewer")
require("gplots")
myCol <- colorRampPalette(c("dodgerblue", "black", "yellow"))(100)
myBreaks <- seq(-3, 3, length.out=101)

par(mar=c(7,4,4,2)+0.1) 
png(filename='Euclidean_Distance.png', width=800, height=750)

#Euclidean distance
heatmap.2(as.matrix(Ribo.filter),
          col=myCol,
          breaks=myBreaks,
          main="Euclidean Distance",
          key=T, keysize=.5,
          scale="none",
          density.info="none",
          reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
          trace="none",
          cexRow=.2,
          cexCol=0.8,
          distfun=function(x) dist(x, method="euclidean"),
          hclustfun=function(x) hclust(x, method="complete"))
graphics.off()


par(mar=c(7,4,4,2)+0.1) 
png(filename='Pearson correlation.png', width=800, height=750)

#1-cor distance
heatmap.2(as.matrix(Ribo.filter),
          col=myCol,
          breaks=myBreaks,
          main="Person Correlation",
          key=T, keysize=1.0,
          scale="none",
          density.info="none",
          reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
          trace="none",
          cexRow=0.2,
          cexCol=0.8,
          distfun=function(x) as.dist(1-cor(t(x))),
          hclustfun=function(x) hclust(x, method="complete"))
graphics.off()


