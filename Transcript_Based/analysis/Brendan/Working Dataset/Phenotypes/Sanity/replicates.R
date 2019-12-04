library(ggpubr)

data("full.pheno.table")
individuals <- full.pheno.table[!duplicated(full.pheno.table$Replicates), "Replicates"]
cor.data <- list()
mean.data <- data.frame(mean = 0, names = 0)
sub.counter <- 0 

# corr between all repticates 
for(i in 1:length(individuals)){
  print(individuals[i])
  sub <- full.pheno.table[which(full.pheno.table$Replicates == individuals[i]),]
  sub.data <- expressed.genes[,which(colnames(expressed.genes) %in% sub[,1])]
  if (length(colnames(sub.data)) == 0){
    print(paste(individuals[i],": skipped. Only has one sample"))
    sub.counter <- sub.counter + 1
  } else{
  cor.table <- cor(na.omit(sub.data), method = "spearman")
  cor.data[[i- sub.counter]] <- cor.table
  names(cor.data)[i- sub.counter] <- individuals[i- sub.counter]
  cor.table[which(cor.table == 1)] <- NA
  mean.data[i- sub.counter,] <- mean(cor.table, na.rm = TRUE)
  mean.data[i- sub.counter,"names"] <- individuals[i- sub.counter]
  print(mean.data)
  }
}


all.cor <- cor(na.omit(expressed.genes), method = "spearman")
all.cor[which(all.cor == 1)] <- NA
mean.data[length(mean.data)+1,] <- mean(all.cor, na.rm = TRUE)
mean.data[which(mean.data$mean == mean(all.cor,na.rm = TRUE)), "names"] <- "All Samples"

ggplot(data = mean.data, aes(x=as.factor(names), y = mean, label = names)) + 
  geom_text()

