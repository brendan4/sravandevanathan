
clusterTranscripts(gene='CHS.771', gown=bg, k=3, method='kmeans')
plotLatentTranscripts(gene='CHS.771', gown=bg, k=3, method='kmeans', returncluster=FALSE)
RPL11 <- collapseTranscripts(gene='CHS.771', gown = bg, meas = "FPKM", method = "kmeans", k = 2)
