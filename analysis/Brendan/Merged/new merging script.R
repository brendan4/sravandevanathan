library(data.table)
library(tidyr)

mergeTables <- function(wd, commonName, colsToMerge){
  #function definations: wd = working directory, commonName = share name of file,
  
  setwd(wd)
  files <- sort(list.files(path ="."))
  firstFile <- files[1]
  
  for(file in files){
    
    setwd(paste(wd,'/',file, sep = ""))
    print(file)
    
    x <- fread(commonName, select = colsToMerge) 
    # will  enter for t_data
    if("gene_name" %in% colnames(x) == TRUE){
      
      names(x)[names(x) == 'gene_name'] <- 'Gene Name'
      names(x)[names(x) == 'start'] <- 'Start'
      names(x)[names(x) == 'end'] <- 'End'
      
      x <- unite(x, `Gene Name`, c(`Gene Name`, Start, End, t_name), remove = TRUE)
      names(x)[names(x) == 'FPKM'] <- paste('FPKM', file)
      
      dups <- c(which(duplicated(x$`Gene Name`)))
      
      #needed for when no dups found 
      if (length(dups) == 0 ){
        xnoDups <- x
      } else {
        xnoDups <- x[-dups] # removes duplicates
      }
    }else {
      x <- unite(x, `Gene Name`, c(`Gene Name`, Start, End), remove = TRUE)
      names(x)[names(x) == 'FPKM'] <- paste('FPKM', file)
      
      xnoDups <- x[-c(which(duplicated(x$`Gene Name`)))] # removes duplicates 
    }
    
    if(file == firstFile){
      assign("mainTable",xnoDups)
      
    }else {
      mainTable<-merge(mainTable, xnoDups,by = "Gene Name", all.x = TRUE, all.y = TRUE)
      assign(file, xnoDups)
    }
  }
  return(mainTable)
}

#gene abundance merge

main.table <- mergeTables(wd = "C:/Users/brendan/Documents/sravandevanathan/ballgown",
                         commonName = "gene_abundance.tab.gz", 
                          colsToMerge = c(2,5,6,8))

main.table <- as.data.frame(main.table) # forms data.table to data.frame
rownames(main.table) <- main.table[,1] #gene namess to row names
main.table <- main.table[,-c(1)] # removes gene names col


boxplot(log2(main.table+.1), 
        names=colnames(main.table), las=2, ylab="log2(FPKM)", 
        main="Distribution of FPKMs for all libraries")

#removal of L2_ACAGTG
main.table <- main.table[,-which(colnames(main.table) %in% c("FPKM L2_ACAGTG"))]

write.table(main.table, "gene_abundance_merged.tab")
GeneAbundance <- read.table("gene_abundance_merged.tab")

# transcript data

main.table <- mergeTables(wd = "C:/Users/brendan/Documents/sravandevanathan/ballgown",
                          commonName = "t_data.ctab.gz", 
                          colsToMerge = c(4,5,6,10,12))

main.table <- as.data.frame(main.table) # forms data.table to data.frame
rownames(main.table) <- main.table[,1] #gene namess to row names
main.table <- main.table[,-c(1)] # removes gene names col


boxplot(log2(main.table+.001), 
        names=colnames(main.table), las=2, ylab="log2(FPKM)", 
        main="Distribution of FPKMs for all libraries")

#removal of L2_ACAGTG
main.table <- main.table[,-which(colnames(main.table) %in% c("FPKM L2_ACAGTG"))]

write.table(main.table,gzfile("transrcipts_merged.ctab.gz")) # writes the table as a .ctab.gz
Transcripts <- read.table(gzfile("transrcipts_merged.ctab.gz")) # will read the table 

