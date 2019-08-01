#'Clean environment 
#'
#'Cleans global enviroment after a long day of working. Always keeps expressed.trans and expressed.genes 
#'
#'@param keep all other objects in environment to keep. Pass objects in as a string (pad with quotes). 
#'Pass multiple objects in as a vector. If object in keep not found in environment a warning will be issued
#'and will terminate the function before removing anything 
#'
#'@return Returns nothing but a clean environment
#'
#'@examples 
#'clean.environment(keep = c("pheno", "expressed.genesN", "gene.results", "bg_filt", "PCA.filt", "filt.set", "filt.names"))
#'clean.environment()
#'clean.enivroment(keep = "bg_filt")

clean.environment <- function(keep = NULL){
  
  if ( length(keep) == 0 ){
    keep <- c("expressed.genes", "expressed.trans", "expressed.genes.GEN", "expressed.trans.GEN")
  }else {
    keep <- c(keep, "expressed.genes", "expressed.trans", "expressed.genes.GEN", "expressed.trans.GEN")
  }
  
  items <- ls(envir=.GlobalEnv)
  
  if (all(keep %in% items) == FALSE){
    not.found = keep[which(!(keep %in% items))]
    warning("Stoping discard- item specified as keep not found. Check: ", not.found)
  }else {
  keep <- c(keep, "discard")
  discard <- items[-c(which(items %in% keep))]
  print(paste("Discarded items: ", discard))
 
  rm(list = c(discard), envir=.GlobalEnv)
  rm(discard)
  }
}

