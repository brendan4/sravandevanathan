
install.packages("tidyr")
library (data.table)
library(tidyr)

#### this section is for t_data only 

setwd("C:/Users/brendan/Documents/sravandevanathan/ballgown")
# set working directory to where files are stored

importing <-c("t_data.ctab")# common name of the files to merge

files <- list.files(path =".")# all the names in the current file path

#merges all the files listed in files aka current path and drops unneed columns
for(file in files){
  setwd(paste("C:/Users/brendan/Documents/sravandevanathan/ballgown/",file, sep = ""))	  
  
  print(file) # which file is being merged
  x <- fread(importing, drop = c(1,2,3,6,7,8,9,11)) # all cols that aren't needed
  x<-unite(x,gene_name,c(gene_name,start,end), remove = TRUE)	# names, start, and end to make unqie name for gene 
  
  names(x)[names(x) == 'FPKM'] <- paste('FPKM', file) # name of the file to the data 
  xnoDups <- x[-c(which(duplicated(x$gene_name)))] # drops all the dups
  
  #indicates the first file where all others will be merged
  if(file == "L2_ACAGTG"){ 
    assign("t_data",xnoDups)
  } else {
    t_data<-merge(t_data, xnoDups,by = "gene_name", all.x = TRUE) 
  }
}

write.table(t_data,gzfile("t_data_merged.ctab.gz")) # writes the table as a .ctab.gz
read.table(gzfile("Merged_t_data.ctab.gz")) # will read the table 

#### this is for gene_aubindance data only 

mergeTables <- function(wd, commonName, firstFile){
  #function definations: wd = working directory, commonName = share name of files, FirstFile = first file to be imported
  
  setwd(wd)
  
  files <- list.files(path =".")
  
  for(file in files){
    setwd(paste(wd,'/',file, sep = ""))
    print(file)
    x <- fread(commonName, drop = c(1,3,4,7,9)) 
    x <- unite(x,`Gene Names`,c(`Gene Name`,Start, End), remove = TRUE)
    names(x)[names(x) == 'FPKM'] <- paste('FPKM', file)
    
    xnoDups <- x[-c(which(duplicated(x$`Gene Names`)))] # removes duplicates 
    
    if(file == firstFile){
      assign("mainTable",xnoDups)
      
    }else {
      mainTable<-merge(mainTable, xnoDups,by = "Gene Names", all.x = TRUE)
      assign(file, xnoDups)
    }
  }
  return(mainTable)
}

mainTable <- mergeTables(wd = "C:/Users/brendan/Documents/sravandevanathan/ballgown",commonName = "gene_abundance.tab", firstFile= "L2_ACAGTG")
# call to function 