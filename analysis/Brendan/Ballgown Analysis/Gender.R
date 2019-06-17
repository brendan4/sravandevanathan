library(dplyr)

pData(bg_filt)=data.frame(id=sampleNames(bg_filt), group=sex)
gene.results=stattest(bg_filt, feature="gene", covariate="Gender", getFC=TRUE, meas="FPKM")
gene.results = data.frame(geneNames=ballgown::geneNames(bg_filt), gene.results)
gene.results = arrange(gene.results,pval)
subset(gene.results,gene.results$qval<0.05)

result_transcripts=stattest(bg_filt, feature="transcript", covariate="Gender", getFC=TRUE, meas="FPKM")
result_transcripts =
  data.frame(geneNames=ballgown::geneNames(bg_filt),
             geneIDs=ballgown::geneIDs(bg_filt), result_transcripts)
result_transcripts = arrange(result_transcripts,pval)
filt.names <- subset(result_transcripts,result_transcripts$qval<0.05)
