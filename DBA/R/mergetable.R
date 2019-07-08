#' Merge Tables 
#' 
#' Merges tables in a directory with a common name 

#' @param wd Working directory in which all folders may be viewed 
#' @param commonName A single string representing the exact name all files share
#' @param colsToMerge integers representing the cols to keep 
#' 
#' @return Returns a new data.table which has combined all the tables into one large table 
#' 
#' @examples 
#' main.table <- mergeTables(wd = "C:/Users/brendan/Documents/sravandevanathan/ballgown",
#' commonName = "gene_abundance.tab", 
#' colsToMerge = c(2,5,6,8))
#' 
#' main.table <- mergeTables(wd = "C:/Users/brendan/Documents/sravandevanathan/ballgown",
#' commonName = "t_data.ctab", 
#' colsToMerge = c(4,5,6,10,12))
#' @import data.table
#' @import tidyr


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
      dups <- c(which(duplicated(x$`Gene Name`)))
      
      if(length(dups) == 0) {
        xnoDups <- x
      }else {
        xnoDups <- x[-dups] # removes duplicates
        }
      xnoDups <- x[-c(which(duplicated(x$`Gene Name`)))] # removes duplicates 
    }
    
    if(file == firstFile){
      assign("mainTable",xnoDups)
      
    }else {
      mainTable <- merge(mainTable, xnoDups,by = "Gene Name", all.x = TRUE, all.y = TRUE)
      assign(file, xnoDups)
    }
  }
  return(mainTable)
}
