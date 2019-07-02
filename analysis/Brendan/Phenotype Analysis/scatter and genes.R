
expressed.genes[grep("EPO", rownames(expressed.genes)),]
test <- expressed.genes[grep("^RBM38",rownames(expressed.genes)),]
plot.var(test)
test <- Var.samples(expressed.genes, gene.list = c("RBM38", "Gapdh") )
test  <- whole.blood[grep("^HB", whole.blood$V1),]


data("full.pheno.table")
diff <- gene.scatter(na.omit(expressed.genes) , 
             "L6_CGATGT", "L6_GGCTAC",
             pheno.table = full.pheno.table,
             names.col = 1, 
             diff.cutoff = 3,
             text.transparency = 1,
             drop.dup.text = TRUE,
             return.gene = TRUE)

test <- expressed.genes[grep("^IL", rownames(expressed.genes)),]

gene.scatter(na.omit(expressed.genes),
             "L6_GCCAAT", "L3_GATCAG", 
             pheno.table = full.pheno.table,
             names.col = 1, 
             diff.cutoff = 3,
             text.transparency = 1,
             drop.dup.text = TRUE,
             return.gene = TRUE)

hemo <- expressed.genes[grep("^HB", rownames(expressed.genes)),]
