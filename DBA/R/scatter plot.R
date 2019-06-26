#' Plot sample vs sample
#' 
#' @import ggplot2
#'
# gene.scatter(expressed.genes, "L6_TTAGGC", "L3_ACTTGA", pheno.table = full.pheno.table, names.col = 1)

gene.scatter <- function (data.set, x.sample, y.sample, 
                          pheno.table = NULL, names.col = NULL, 
                          text.transparency = .5, point.transparency = .1,
                          min.cutoff = log(.101), diff.cutoff = 2,
                          drop.dup.text = FALSE){
  #data prep
  data.set <- log(na.omit(data.set)+0.1)
  X = data.set[,x.sample]
  Y = data.set[,y.sample]
  data.set <- data.set[,which(colnames(data.set) 
                              %in% c(x.sample,y.sample))]
  #phenotype info
  if(is.null(pheno.table) == FALSE){
    if(is.null(names.col) == TRUE){
      warning("must provide which column holds the sample names in phenotype table as names.col")
    } else{
      sub <- pheno.table[which(pheno.table[,names.col] 
                               %in% c(x.sample, y.sample)),]
      print(sub)
    }
  }
  #subseting data for teex labels
  subdata <- subset(data.set, ((data.set[,x.sample] - data.set[,y.sample]) > diff.cutoff 
                               | (data.set[,x.sample] - data.set[,y.sample]) < -diff.cutoff) 
                    & !(data.set[y.sample] < min.cutoff | data.set[x.sample] < min.cutoff))
  #pretty gene names 
  if(drop.dup.text == TRUE){
    subdata <- pretty.gene.name(subdata, as.row.names = TRUE, remove.dups = TRUE)
  }else {
  subdata <- pretty.gene.name(subdata)
    if(sum(duplicated(subdata$pretty)) == 0){
      rownames(subdata) <- subdata$pretty
    }else{
      warning("Pretty gene name created duplicates: using long gene names instead")
    }
  }
  
  #ggplots data
  ggplot(data.set, aes(x = X, y = Y)) + 
    geom_point(shape = 1, alpha = point.transparency) +
    theme_minimal() +
    xlab(x.sample) +
    ylab(y.sample) +
    scale_colour_hue(l = 50) + # Use a slightly darker palette than normal
    geom_smooth(method = lm,   # Add linear regression lines
                se = TRUE,    # Don't add shaded confidence region
                fullrange=TRUE)+ # Extend regression lines
    geom_point(shape = 23,
               fill = "red", 
               data = subdata, 
               aes(x = subdata[,x.sample],
                   y = subdata[,y.sample]),
               alpha = .5) +
    geom_text(data = subdata, 
              aes(x = subdata[,x.sample],
                  y = subdata[,y.sample], 
                  label = row.names(subdata)),
              show.legend = FALSE,
              position = "jitter", 
              alpha = text.transparency)
              
}


#gene.scatter(expressed.genes, "L6_TTAGGC", "L3_ACTTGA", pheno.table = full.pheno.table, names.col = 1)

#subdata <- data.set
#subdata$var <- subdata[,x.sample] - subdata[,y.sample]
#subset(subdata, subdata$var > 2)

# subdata <- subset(data.set, (5 > data.set[,x.sample] & data.set[,y.sample] > 6) | (5 > data.set[,y.sample] & data.set[,x.sample] > 6))