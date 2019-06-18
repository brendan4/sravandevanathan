RPL11 <- gene.info[grep("^RPL11",gene.info[,2]),] # gene_id from plotTransript

plotTranscripts(gene='CHS.771', gown=bg, samples='L2_GCCAAT', 
                meas='FPKM', colorby='transcript', 
                main='transcripts from gene XLOC_000454: L2_GCCAAT, FPKM')

all.plotTranscripts <- function(Gene, samples.list = NULL, gene.info, bg){
  pattern = gene.info[grep(paste("^", Gene, sep = ""),gene.info$gene_name),]
  #if (pattern == 0 ){
  # print(paste(pattern,": failed to find match", sep = ""))
  # break else {
  geneID <- as.character(pattern[1,1]) # uses first in data.frame as geneID
  #if (samples.list == NULL) {
  #samples.list <-  sampleNames(bg)
  #}
  for (Sample in samples.list){
    title <-  paste("transcripts from",Gene, ":", Sample, "FPKM", sep = " ")
    file <- paste(Sample," - ",Gene,".png", sep ="")
    par(mar=c(7,4,4,2)+0.1) 
    png(filename = file, width=650, height=450)
    plotTranscripts(gene= geneID, gown=bg, samples= Sample, 
                    meas='FPKM', colorby='transcript', 
                    main=  title)
    graphics.off()
    
  }
}


setwd("C:\\Users\\brendan\\Documents\\sravandevanathan\\analysis\\Brendan\\Ballgown_analysis")
all.plotTranscripts(Gene = "RPL11", samples.list = sampleNames(bg), gene.info = gene.info, bg = bg )
