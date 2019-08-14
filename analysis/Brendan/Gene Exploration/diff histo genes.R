data("expressed.genes.GEN")
data("pheno.more.samples")
CDK <- t(expressed.genes.GEN[grep("CDK11A", rownames(expressed.genes.GEN)), ])
pheno.more.samples$CDK <- CDK


summary(pheno.more.samples$CDK[which(pheno.more.samples$pheno == "C")])
c <- pheno.more.samples[which(pheno.more.samples$pheno == "C"), ]
out <- c[which(c$CDK < 32),]
ggplot(pheno.more.samples, aes(x = pheno, y= CDK, color = pheno))+
  geom_boxplot()+
  labs(title="CDK11A FPKM Expression",x="Phenotype", y = "CDK11A FPKM")+
  geom_text(data = out, aes(x = pheno, y = CDK, label =`colnames(expressed.genes)`))


pheno.more.samples[which(pheno.more.samples$pheno == "SS"), "pheno"] <- "S"
pheno.more.samples<- gene.histo("RPL11",
                                data.set = expressed.genes.GEN, 
                                pheno.data = pheno.more.samples)

ggplot(pheno.more.samples, aes(x = pheno, y= gene, color = pheno))+
  geom_boxplot()+
  labs(title="RPL11 FPKM Expression",x="Phenotype", y = "RPL11 FPKM")

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6416817/
DBA.related.Ribo <- c("RPL5", "RPL11", "RPL35A", "RPS7", "RPS10", "RPS17", "RPS19", "RPS24", "RPS26",
                      "RPL3", "RPL7", "RPL9", "RPL14", "RPL19", "RPL23A", "RPL26", "RPL35", "RPL36", "RPS8"
                      ,"RPS15", "RPS27A", "RPL18","RPL35")

DBA <- filter.genes(expressed.genes.GEN, DBA.related.Ribo)
