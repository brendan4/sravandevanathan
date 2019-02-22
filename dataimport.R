install.packages("tidyr")
library (data.table)
library(tidyr)
setwd("C:/Users/brendan/Documents/sravandevanathan/ballgown")

importing <-c("gene_abundance.tab")

files <- list.files(path =".")


for(file in files){
  setwd(paste("C:/Users/brendan/Documents/sravandevanathan/ballgown/",file, sep = ""))
  print(file)
  x <- fread(importing, drop = c(1,3,4,7,9)) 
  x <- unite(x,`Gene Names`,c(`Gene Name`,Start, End), remove = TRUE)
  names(x)[names(x) == 'FPKM'] <- paste('FPKM', file)
  
  if(file == "L2_ACAGTG"){
    assign("mainTable",x)
    } else {
    mainTable<-merge(mainTable, x,by = "Gene Names", all.x = TRUE)
    #mainTable <- na.omit(mainTable)
    assign(file, x)
  }
}


merge_data <- function(mainTable, merger){
  merge(mainTable, merger, all.x = TRUE)
} 
merge_data(mainTable,L2_ACTTTGA)
mainTable<-merge(mainTable, L2_ACTTTGA, by.x = "Gene Names", by.y = "Gene Names",all.x = TRUE, suffixes = files[1])
