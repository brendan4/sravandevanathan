#' Plot sample vs sample
#' 
#' @example gene.scatter(expressed.genes, "L6_TTAGGC", "L3_ACTTGA", pheno.table = full.pheno.table, names.col = 1)

gene.scatter <- function (data.set, x.sample, y.sample, pheno.table = NULL, names.col = NULL){
  data.set <- na.omit(data.set)
  x = data.set[,x.sample]
  y = data.set[,y.sample]
  
  
  if(is.null(pheno.table) == FALSE){
    if(is.null(names.col) == TRUE){
      warning("must provide which column holds the sample names in phenotype table as names.col")
    } else{
      sub <- pheno.table[which(pheno.table[,names.col] %in% c(x.sample, y.sample)),]
      print(sub)
    }
  }
  ggplot(expressed.genesNA, aes(x=x, y=y)) + geom_point(shape=1) +
    scale_colour_hue(l=50) + # Use a slightly darker palette than normal
    geom_smooth(method=lm,   # Add linear regression lines
                se=FALSE,    # Don't add shaded confidence region
                fullrange=TRUE) # Extend regression lines
  
}



