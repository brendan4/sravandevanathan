library(dplyr)
library(propr)
library(circlize)
library(ComplexHeatmap) 
#raw data 
setwd("~/sravandevanathan/analysis/Brendan/Working Dataset/New Annotations")
GeneAbundance <- read.table("gene_abundance_merged.tab")
Transcripts <- read.table(gzfile("transrcipts_merged.ctab.gz"))

#https://www.genenames.org/data/genegroup/#!/group/728s
s.ribo <- read.csv("C:/Users/brendan/Documents/sravandevanathan/analysis/Brendan/Gene Exploration/sribo.csv",
                   header = T, stringsAsFactors = F, skip = 1)
l.ribo <- read.csv("C:/Users/brendan/Documents/sravandevanathan/analysis/Brendan/Gene Exploration/lribo.csv",
                   header = T, stringsAsFactors = F, skip = 1)
all.ribo <- c(s.ribo$Approved.symbol, l.ribo$Approved.symbol)

#data splitting
data("pheno.more.samples")
carr <- pheno.more.samples %>% filter(pheno == "C")
wild <- pheno.more.samples %>% filter(pheno == "W")
ser <- pheno.more.samples %>% filter(pheno == "S"| pheno == "SS")
carr[duplicated(carr$Replicates),]

#sanity options
sanity <- carr[c(5,7,4,9),]
sanity <- carr[c(4,10,9,5),]
sanity <- pheno.more.samples[c(5,6,7,12, 27, 31, 1,4,3,24),]
sanity <- pheno.more.samples %>% filter(pheno == "S")

#options for sliceing
carr.genes <- GeneAbundance[,which(colnames(expressed.genes.GEN) %in% carr$`colnames(expressed.genes)`)]
wild.genes <- GeneAbundance[which(colnames(expressed.genes.GEN) %in% wild$`colnames(expressed.genes)`)]
ser.genes <- GeneAbundance[which(colnames(expressed.genes.GEN) %in% ser$`colnames(expressed.genes)`)]
sanity.genes <- GeneAbundance[which(colnames(GeneAbundance) %in% sanity$`colnames(expressed.genes)`)]

# fin tune the sum function for filtering down the line 
all.genes <- c(all.ribo,"ATF1", "MYC", "TP53", "CDK11a", "HBB", 'NOB1', "TCTA") #housekeeping,
pheno.filt <- filter.genes(sanity.genes, all.genes, lazy = F)
pheno.filt <- pheno.filt[,-which(colnames(pheno.filt) %in% "pretty")]
pheno.filt<- t(pheno.filt)
pheno.more.samples[which(pheno.more.samples$`colnames(expressed.genes)` %in% rownames(pheno.filt)),] #sanity check
keep.test <- apply(pheno.filt, 2, function(x) sum(x >= 20) >= 4)
keep.test[which(names(keep.test) %in% 
             c("ATF1_50764101_50821162", "MYC_127735434_127742951", "TP53_7661779_7687550", "CDK11A_1702379_1724357", 'TCTA_49412206_49416475', 'NOB1_69741871_69754926'))] <- TRUE
keep.test[which(names(keep.test) %in% c("RPS4Y1_2841602_2867268"))] <- F
ribo.sub <- names(keep.test)[which(keep.test == T)]


#unnorm data and feed to propr
unnorm.expressed.genes <- unnormilize(carr.genes, reads)
unnorm.expressed.genes <- t(unnorm.expressed.genes)

#making cutoffs for select
keep <- apply(unnorm.expressed.genes, 2, function(x) sum(x >= 75) >= 10) 
#marking all ribo related as true
keep[which(names(keep) %in% ribo.sub)] <- TRUE
rho <- propr(unnorm.expressed.genes, metric = "rho", select = keep) # pror functioncall

#rho data processing
rho.DBA <- subset(rho, select = ribo.sub)
DBA.matrix <- rho.DBA@matrix
rownames(DBA.matrix) <- str_extract(rownames(DBA.matrix), "^[-a-zA-Z0-9\\.]*")
colnames(DBA.matrix) <- str_extract(colnames(DBA.matrix), "^[-a-zA-Z0-9\\.]*")

#heatmap ploting
Heatmap(DBA.matrix,
        column_names_side = "bottom",
        column_names_gp = gpar(fontsize = 8),
        row_names_gp = gpar(fontsize = 8),
        row_names_side = "left")#col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))

prism(rho.DBA, plotly = T, k=3)
clean.environment(keep = c("unnormilize", "reads"))


keep <- colnames(expressed.genes.GEN) %in% c("L6_CGATGT", 'L2_CTTGTA', "L2_GCCAAT", 'L3_CGATGT')
counts <- expressed.genes.GEN[,keep]
counts <- t(counts)
counts <- counts[,1:100]
group <- ifelse(rownames(counts) %in% c("L6_CGATGT", 'L2_CTTGTA'), "S", "C")
pd <- propd(counts, group, alpha = NA, p = 10)
look <- pd@results

pd.d <- setDisjointed(pd)
tab <- getResults(pd.d)
plot(pd.d@counts[, 39], pd.d@counts[, 37], col = ifelse(pd.d@group == "S", "red", "blue"))
grp1 <- pd.d@group == "S"
grp2 <- pd.d@group != "S"
abline(a = 0, b = pd.d@counts[grp1, 37] / pd.d@counts[grp1, 39], col = "red")
abline(a = 0, b = pd.d@counts[grp2, 37] / pd.d@counts[grp2, 39], col = "blue")
plot(pd.d@counts[, 37] / pd.d@counts[, 39],
     col = ifelse(pd.d@group == "S", "red", "blue"))