library(propr)

expressed.genes2 <- na.omit(expressed.genes)
expressed.genes2 <- pretty.gene.name(expressed.genes2, as.row.names = T, remove.dups = T)
expressed.genes2 <- t(expressed.genes2)
keep <- apply(expressed.genes2, 2, function(x) sum(x >= 10) >= 10)

# data prep
hemo <- filter.genes(expressed.genes, c("HBB", "HBG1", "HBG2", "HBA2", "HBA1") )
hemo <- pretty.gene.name(hemo, as.row.names = TRUE)
hemo <- t(hemo)

phi <- propr(expressed.genes2, metric = "rho", select = keep)
updateCutoffs(phi, cutoff =  seq(.05, .95, .3))

#Subset option 1
phiHBB.HBG1 <- subset(phi, select = c("HBB", "HBG1"))
phiHBB.HBG1@matrix

#subset option 2
phi99 <- phi[">", .95]
simplify(phi99)
phi99@matrix

#ploting options 
plot(phi99)
snapshot(phi, prompt = TRUE, plotly = TRUE)
pca(phi, group = pheno$pheno , prompt = TRUE, plotly = FALSE) #group
dendrogram(phi, prompt = TRUE, plotly = TRUE)
smear(phi99, prompt = TRUE, plotly = TRUE)

#Rho plots 
slate(phi, k, prompt = TRUE, plotly = FALSE)
bucket(phi, group, k = 5, prompt = TRUE, plotly = FALSE)
prism(phi, k = 5, prompt = TRUE, plotly = FALSE)
bokeh(phi, k, prompt = TRUE, plotly = FALSE)
cytescape(phi99, col1, col2, prompt = TRUE, d3 = FALSE)


