diff.genes.cleanup <- function(diff.results, gown, subset = FALSE){
  
  #gene names matching
  indices <- match(diff.results$id, texpr(gown, 'all')$gene_id)
  gene_names_for_result <- texpr(gown, 'all')$gene_name[indices]
  diff.results <- data.frame(geneNames = gene_names_for_result, diff.results)
  
  if(subset == TRUE){

  #processing 
  diff.results = arrange(diff.results, pval)
  filt.names <- subset(diff.results, diff.results$qval<0.05)
  
  #filters outs periods
  filt.set <- filter.out.genes(filt.names, gene.list = c("\\."), by.rownames = FALSE, col = 1)
  diff.genes <- filter.genes(expressed.genes, filt.set$geneNames)
  return(diff.genes)
  }
  return(diff.results)
}