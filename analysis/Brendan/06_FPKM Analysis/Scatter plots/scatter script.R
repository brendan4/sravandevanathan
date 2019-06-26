 
#L6_CGATGT -> IV.5  vs  L6_GGCTAC -> III.13 (some symptoms)
gene.scatter(expressed.genes, "L6_CGATGT", "L6_GGCTAC", 
             pheno.table = full.pheno.table, names.col = 1,
             diff.cutoff = 3,
             text.transparency = .75, 
             point.transparency = .35)

#L6_CGATGT -> IV.5  vs  L3_CAGATC III.11 (no mutation)
gene.scatter(expressed.genes, "L6_CGATGT", "L3_CAGATC", 
             pheno.table = full.pheno.table, names.col = 1,
             diff.cutoff = 3,
             min.cutoff = log(.3),
             text.transparency = .75, 
             point.transparency = .35)

# IV.5 replicate (L6_CGATGT  vs  L2_CTTGTA)
gene.scatter(expressed.genes, "L6_CGATGT", "L2_CTTGTA", 
             pheno.table = full.pheno.table, names.col = 1,
             min.cutoff = log(1),
             diff.cutoff = 2, 
             text.transparency = 1, 
             point.transparency = .35)
cor(na.omit(expressed.genes[,c("L6_CGATGT", "L2_CTTGTA")]),method = "spearman")

# I.1 repliacte (L6_TTAGGC  vs  L3_ACTTGA) - maybe not the same person 
gene.scatter(expressed.genes, "L6_TTAGGC", "L3_ACTTGA", 
             pheno.table = full.pheno.table, names.col = 1,
             min.cutoff = log(1),
             diff.cutoff = 2, 
             text.transparency = .8, 
             point.transparency = .35)
cor(na.omit(expressed.genes[,c("L6_TTAGGC", "L3_ACTTGA","L6_CGATGT", "L2_CTTGTA")]), method = "spearman")
