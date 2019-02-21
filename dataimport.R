library (data.table)
setwd("C:/Users/brendan/Documents/sravandevanathan/ballgown")

importing <-c("gene_abundance.tab")
files <- list.files(path =".")
  
for(file in files){
  setwd(paste("C:/Users/brendan/Documents/sravandevanathan/ballgown/",file, sep = ""))
  x <- fread(importing, drop = c(1,3,4,7,9)) 
  assign(file, x)
}

for(table in files){
  
}


testmerge <- merge(test1, test2, by.x = c(1,2), by.y = c(1,2), all.x = TRUE)
