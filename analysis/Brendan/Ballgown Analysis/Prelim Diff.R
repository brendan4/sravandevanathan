library(dplyr)

pData(bg_filt) = data.frame(id = sampleNames(bg_filt), group = pheno)
gene.results = stattest(bg_filt, feature="gene", covariate="group.pheno", getFC=TRUE, meas="FPKM")

#gene names 
indices <- match(gene.results$id, texpr(bg_filt, 'all')$gene_id)
gene_names_for_result <- texpr(bg_filt, 'all')$gene_name[indices]
gene.results <- data.frame(geneNames = gene_names_for_result, gene.results)

#processing 
gene.results = arrange(gene.results, pval)
filt.names <- subset(gene.results, gene.results$qval<0.05)

#filter out . from data set as it acts as a regular expression 
filt.set <- filter.out.genes(filt.names, gene.list = c("\\."), by.rownames = FALSE, col = 1)
diff.genes <- filter.genes(expressed.genes, filt.set$geneNames)

