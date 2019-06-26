

gene.scatter <- function (data.set, x.sample, y.sample, 
                          pheno.table = NULL, names.col = NULL, 
                          text.transparency = .5, point.transparency = .1,
                          min.cutoff = log(.101), diff.cutoff = 2){
  data.set <- log10(na.omit(data.set)+0.1)
  X = data.set[,x.sample]
  Y = data.set[,y.sample]
  data.set <- data.set[,which(colnames(data.set) 
                              %in% c(x.sample,y.sample))]
  
  
  if(is.null(pheno.table) == FALSE){
    if(is.null(names.col) == TRUE){
      warning("must provide which column holds the sample names in phenotype table as names.col")
    } else{
      sub <- pheno.table[which(pheno.table[,names.col] 
                               %in% c(x.sample, y.sample)),]
      print(sub)
    }
  }
  subdata <- subset(data.set, ((data.set[,x.sample] - data.set[,y.sample]) > diff.cutoff 
                               | (data.set[,x.sample] - data.set[,y.sample]) < -diff.cutoff) 
                    & !(data.set[y.sample] < min.cutoff | data.set[x.sample] < min.cutoff))
  
  ggplot(data.set, aes(x = X, y = Y)) + 
    geom_point(shape=1, alpha = point.transparency) +
    theme_minimal() +
    xlab(x.sample) +
    ylab(y.sample) +
    scale_colour_hue(l=50) + # Use a slightly darker palette than normal
    geom_smooth(method=lm,   # Add linear regression lines
                se=TRUE,    # Don't add shaded confidence region
                fullrange=TRUE)+ # Extend regression lines
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