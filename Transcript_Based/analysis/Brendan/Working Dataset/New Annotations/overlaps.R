clean.environment()
overlap <- intersect(rownames(expressed.genes.GEN), rownames(expressed.genes.more))
non.overlap.GEN <- expressed.genes.GEN[which(!(rownames(expressed.genes.GEN) %in% overlap)),]
non.overlap.old <- expressed.genes.more[which(!(rownames(expressed.genes.more) %in% overlap)),]
non.overlap.GEN <- pretty.gene.name(non.overlap.GEN)
non.overlap.old <- pretty.gene.name(non.overlap.old)
non.overlap.GEN <- non.overlap.GEN[!duplicated(non.overlap.GEN$pretty),]
non.overlap.old <- non.overlap.old[!duplicated(non.overlap.old$pretty),]
non.overlap.GEN <- non.overlap.GEN[which(non.overlap.GEN$pretty %in% non.overlap.old$pretty),]
non.overlap.old <- non.overlap.old[which(non.overlap.old$pretty %in% non.overlap.GEN$pretty),]
comp <- data.frame(new = rownames(non.overlap.GEN), old = rownames(non.overlap.old))


non.overlap.GEN <- pretty.gene.name(non.overlap.GEN, as.row.names = TRUE, remove.dups = TRUE)
non.overlap.old <- pretty.gene.name(non.overlap.old, as.row.names = TRUE, remove.dups = TRUE)
pretty <- intersect(rownames(non.overlap.GEN), rownames(non.overlap.old))
          