
gene.histo <- function(gene, data.set, pheno.data){
  gene.data <- data.set[grep(paste("^", gene, sep = ""), rownames(data.set)), ]
  rownames(gene.data) <- str_extract(rownames(gene.data), "^[-a-zA-Z0-9\\.]*")
  gene.data <- gene.data[grep(paste("^", gene,"$", sep = ""), rownames(gene.data)),]
  gene.data <- t(gene.data)
  pheno.data["gene"] <- gene.data
  
  return(pheno.data)
}

