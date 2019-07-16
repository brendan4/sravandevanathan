library(propr)

#for clr
expressed.genes2 <- na.omit(expressed.genes)
expressed.genes2 <- pretty.gene.name(expressed.genes2, as.row.names = T, remove.dups = T)
expressed.genes2 <- t(expressed.genes2)
keep <- apply(expressed.genes2, 2, function(x) sum(x >= 10) >= 10)

# data prep for alr 
hemo <- c("HBB", "HBG1", "HBG2", "HBA2", "HBA1", "HSPA4", "PPID") 
hemo <- filter.genes(expressed.genes, hemo )
hemo <- pretty.gene.name(hemo, as.row.names = TRUE)
hemo <- t(hemo)

phi <- propr(expressed.genes2, metric = "rho", select = keep)
updateCutoffs(phi, cutoff =  seq(.05, .95, .3))

phi <- propr(hemo, metric = "rho", ivar = 8)
updateCutoffs(phi, cutoff =  seq(.05, .95, .3))

#Subset option 1
phiHBB.HBG1 <- subset(phi, select = hemo)
phiHBB.HBG1@matrix

#subset option 2
phi99 <- phi[">", .95]
simplify(phi99)
phi99@matrix

#ploting options 
plot(phiHBB.HBG1)
snapshot(phiHBB.HBG1, prompt = TRUE, plotly = TRUE)
pca(phiHBB.HBG1, group = pheno$pheno , prompt = TRUE, plotly = TRUE) #group
dendrogram(phiHBB.HBG1, prompt = TRUE, plotly = TRUE)
smear(phiHBB.HBG1, prompt = TRUE, plotly = TRUE)

#Rho plots 
slate(phi, k, prompt = TRUE, plotly = FALSE)
bucket(phi99, group = pheno$pheno , k = 5, prompt = TRUE, plotly = FALSE)
prism(phi99, k = 5, prompt = TRUE, plotly = TRUE)
bokeh(phi, k = 5, prompt = TRUE, plotly = FALSE)
cytescape(phi99, prompt = TRUE)


smear(phi99, prompt = T, plotly = T)
pca(phi99, group = pheno$pheno)
snapshot(phi99, plotly = T)
