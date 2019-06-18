data("pheno.basic")
data("sex")

pData(bg_filt) = data.frame(id = sampleNames(bg_filt), group = pheno.basic$pheno, group = sex$Gender)
gene.results = stattest(bg_filt, feature="gene", covariate="group", getFC=TRUE, meas="FPKM", adjustvars = "group.1")

diff.genes <- diff.genes.cleanup(gene.results,bg_filt,subset = TRUE)
