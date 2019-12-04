validation.helper <- function(full.pheno, data.set, class, cate, sub = NULL){
  counter <- 0
  mutant.flag = FALSE
  if(is.null(sub) == FALSE) {mutant.flag = TRUE}
 
  if(mutant.flag == FALSE){
    id.holder <- unlist(full.pheno[[class]][[cate]], 
                        use.names = FALSE)
  } else{
    id.holder <- unlist(full.pheno[[class]][[cate]][sub], 
                        use.names = FALSE)
  }
  
  counter <- length(id.holder)
  lost<- which(!(id.holder %in% colnames(expressed.genes)))
  
  if(length(lost) == 0){
    num <- 0
    lost <- 0
    if(mutant.flag == FALSE){
      print(paste(names(full.pheno[[class]][cate]), 
                  ":", full.pheno[[class]][cate], ": is good"), sep= " ")
    } else{
      print(paste(names(full.pheno[[class]][[cate]][sub]),
                  ":", full.pheno[[class]][[cate]][sub], ": is good"), sep= " ")
    }
  } else{
    num <- length(lost)
    if(mutant.flag == FALSE){
      lost <- full.pheno[[class]][[cate]][c(lost)]
    } else{
      lost <- full.pheno[[class]][[cate]][[sub]][c(lost)]
    }
  }
  items <- list(num, c(lost), counter)
  return(items)
}

validation <- function(full.pheno, data.set){
  num.lost <- 0
  lost.values <- c()
  counter <- 0
  
  for(class in 1:length(full.pheno)){
    
    if (names(full.pheno)[class] == "Wildtype"){
      print("Checking wildtype data")
      for(cate in 1:length(full.pheno[[class]])){
  
        items <- validation.helper(full.pheno, data.set, class, cate)
        num <- items[[1]]
        counter <- counter + items[[3]]
      
        if(num == 0){
          
        } else{
          lost <- items[[2]]
          lost.values <- append(lost.values, lost)
          num.lost <- num.lost + num
        }
      }
    } else{
      print("Checking Mutant data: ")
      for(cate in 1:length(full.pheno[[class]])){
        for(sub in 1:length(full.pheno[[class]][[cate]])){
          items <- validation.helper(full.pheno, data.set, class, cate, sub = sub)
          num <- items[[1]]
          counter <- counter + items[[3]]
          
          if(num == 0){
            
          } else{
            lost <- items[[2]]
            lost.values <- append(lost.values, lost)
            num.lost <- num.lost + num
          }
        }
        
      }
    }
  }
  print("OVer all results: ")
  print(paste("Number of samples lost: ",num.lost))
  print("Sample tags: ")
  print(lost.values)
  print(paste("Out of a total: ", counter))
  return(lost.values)
}

lost.samples<- validation(full.pheno, expressed.genes)  

#L3_ATCACG thrown out 
lost.samples[which(lost.samples %in% colnames(expressed.genesN))]
lost.samples <- lost.samples[-which(lost.samples %in% colnames(expressed.genesN))]

#L2_ACAGTG found in file listing: likely removed early on: COMFIRmED BY LOOKING AT MERGED
setwd("C:/Users/brendan/Documents/sravandevanathan/Ballgown")
length(list.files()) - 3
lost.samples[which(lost.samples %in% list.files())]

#nothing found 
library(data.table)
setwd("C:/Users/brendan/Documents/sravandevanathan/analysis/Brendan/Working Dataset/01_Raw Data")
raw_data <- fread("gene_abundance_merged.tab")

lost.samples[which(paste("FPKM",lost.samples) %in% colnames(raw_data))]


#L6_ACAGTC remains unaccounted for 
##UPDATE: L6_ACAGTG is file - it was mispelled