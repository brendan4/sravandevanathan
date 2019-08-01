#processed data sanity check

data(expressed.genes.GEN)
data("expressed.trans.GEN")
data("expressed.trans")
data("expressed.genes.more")
data("expressed.genes.GEN")
data("pheno.more.samples")



expressed.genes.GEN <- expressed.genes.GEN[,c(5:36,1,2,3,4)]
expressed.genes.more <- na.omit(expressed.genes.more)
colnames(expressed.genes.GEN)[which(colnames(expressed.genes.GEN) 
                                    %in% c("LIB1.98408313", "LIB3.98406315", "LIB5.98403313", "LIB6.98398316"))] <- c("LIB1-98408313","LIB3-98406315","LIB5-98403313","LIB6-98398316")

#checking colnames intersection
intersect(colnames(expressed.genes.more), colnames(expressed.genes.GEN))
tes <- colnames(expressed.genes.GEN)[which(colnames(expressed.genes.GEN) %in% colnames(expressed.genes.more))]

#gene names intersections
expressed.genes.more <- pretty.gene.name(expressed.genes.more)
expressed.genes.GEN <- pretty.gene.name(expressed.genes.GEN)

#intersecting data frames
t3 <- intersect(expressed.genes.more$pretty, expressed.genes.GEN$pretty)
inta.GEN <- expressed.genes.GEN[which(expressed.genes.GEN$pretty %in% t3), ]
inta.more <- expressed.genes.more[which(expressed.genes.more$pretty %in% t3), ]

#a look at dups
dups<- inta.GEN[duplicated(inta.GEN$pretty),"pretty"]
dups.extended <- inta.GEN[which(inta.GEN$pretty %in% dups),]
dups<- inta.more[duplicated(inta.more$pretty),"pretty"]
dups.extended <- inta.more[which(inta.more$pretty %in% dups),]

#handling non-overlapping gene regions 
inta.more <- inta.more[, -which(colnames(inta.more) %in% "pretty") ]
inta.GEN <- inta.GEN[, -which(colnames(inta.GEN) %in% "pretty") ]

inta.more <- pretty.gene.name(inta.more, as.row.names = TRUE, remove.dups = TRUE)
inta.GEN <- pretty.gene.name(inta.GEN, as.row.names= TRUE, remove.dups = TRUE)
cor.table <- cor(inta.more, inta.GEN, method = "spearman")

clean.environment()

