install.packages("tidyr")
library (data.table)
library(tidyr)
setwd("C:/Users/brendan/Documents/sravandevanathan/ballgown")

importing <-c("gene_abundance.tab")
files <- list.files(path =".")


for(file in files){
  setwd(paste("C:/Users/brendan/Documents/sravandevanathan/ballgown/",file, sep = ""))
  x <- fread(importing, drop = c(1,3,4,7,9)) 
  x<-unite(x,`Gene Names`,c(`Gene Name`,Start,End), remove = TRUE)
  
  assign(file, x)
  
}


