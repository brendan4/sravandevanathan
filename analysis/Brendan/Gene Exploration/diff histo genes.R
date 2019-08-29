data("expressed.genes.GEN")
data("pheno.more.samples")

#CDK11a histo plot
CDK <- t(expressed.genes.GEN[grep("CDK11A", rownames(expressed.genes.GEN)), ])
pheno.more.samples$CDK <- CDK

summary(pheno.more.samples$CDK[which(pheno.more.samples$pheno == "C")])
c <- pheno.more.samples[which(pheno.more.samples$pheno == "C"), ]
out <- c[which(c$CDK < 32),]
ggplot(pheno.more.samples, aes(x = pheno, y= CDK, color = pheno))+
  geom_boxplot()+
  labs(title="CDK11A FPKM Expression",x="Phenotype", y = "CDK11A FPKM")+
  geom_text(data = out, aes(x = pheno, y = CDK, label =`colnames(expressed.genes)`))

#histo for RPL11 
#forming groups and data for gene
pheno.more.samples[which(pheno.more.samples$pheno == "SS"), "pheno"] <- "S"
pheno.more.samples<- gene.histo("RPL11",
                                data.set = expressed.genes.GEN, 
                                pheno.data = pheno.more.samples)

group.means <- data.frame()
group.means["C","mean"] <- mean(pheno.more.samples[which(pheno.more.samples$pheno == "C"), "gene"])
group.means["S","mean"] <- mean(pheno.more.samples[which(pheno.more.samples$pheno == "S"), "gene"])
group.means["W","mean"] <- mean(pheno.more.samples[which(pheno.more.samples$pheno == "W"), "gene"])
group.means["pheno"] <- rownames(group.means)

ggplot(pheno.more.samples, aes(x = pheno, y = gene, color = pheno))+
  geom_point()+
  geom_bar(data = group.means, aes(x= pheno, y = mean), stat= "identity", alpha = .2)+
  labs(title="RPL11 FPKM Expression",x="Phenotype", y = "RPL11 FPKM")+
  theme_classic()

#DBA releated ribo genes
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6416817/
DBA.related.Ribo <- c("RPL5", "RPL11", "RPL35A", "RPS7", "RPS10", "RPS17", "RPS19", "RPS24", "RPS26",
                      "RPL3", "RPL7", "RPL9", "RPL14", "RPL19", "RPL23A", "RPL26", "RPL35", "RPL36", "RPS8"
                      ,"RPS15", "RPS27A", "RPL18","RPL35")

DBA <- filter.genes(expressed.genes.GEN, DBA.related.Ribo, lazy = F)

plot.FPKM <- function(data.set, pheno){
  
  for(i in 1:nrow(data.set)){
    print(rownames(data.set)[i])
    pheno.more.samples <- gene.histo(rownames(data.set)[i],
                                    data.set = data.set, 
                                    pheno.data = pheno)
    print(rownames(data.set))
    group.means <- data.frame()
    group.means["C","mean"] <- mean(pheno.more.samples[which(pheno.more.samples$pheno == "C"), "gene"])
    group.means["S","mean"] <- mean(pheno.more.samples[which(pheno.more.samples$pheno == "S"), "gene"])
    group.means["W","mean"] <- mean(pheno.more.samples[which(pheno.more.samples$pheno == "W"), "gene"])
    group.means["pheno"] <- rownames(group.means)
    
    ggplot(pheno.more.samples, aes(x = pheno, y = gene, color = pheno))+
      geom_point()+
      geom_bar(data = group.means, aes(x= pheno, y = mean), stat= "identity", alpha = .2)+
      labs(title="RPL11 FPKM Expression",x="Phenotype", y = "RPL11 FPKM")+
      theme_classic()
    pheno.more.samples <- pheno.more.samples[,-which(colnames(pheno.more.samples) %in% "gene")]
  }
}

plot.FPKM(DBA, pheno.more.samples)

