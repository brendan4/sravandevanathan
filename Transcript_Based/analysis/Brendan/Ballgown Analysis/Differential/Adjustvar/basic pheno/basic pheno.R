data("pheno.basic")
data("sex")

pData(bg_filt) = data.frame(id = sampleNames(bg_filt), 
                            group = pheno.basic$pheno)

gene.results = stattest(bg_filt, 
                        feature="gene", 
                        covariate="group", 
                        getFC=TRUE, 
                        meas="FPKM")

#clean up data
diff.genes <- diff.genes.cleanup(gene.results,
                                 bg_filt, 
                                 subset = TRUE)
