clean.environment <- function( keep = NULL ){
  if ( length(keep) == 0 ){
    keep <- c("expressed.genes", "expressed.trans","discard","clean.environment")
  }else {
    keep <- c(keep, "expressed.genes", "expressed.trans","discard","clean.environment")
  }
  
  items <- ls(envir=.GlobalEnv)
  print(items)
  discard <- items[-c(which(items %in% keep))]
  print(discard)
 
  rm(list = c(discard), envir=.GlobalEnv)
  rm(discard)
}

