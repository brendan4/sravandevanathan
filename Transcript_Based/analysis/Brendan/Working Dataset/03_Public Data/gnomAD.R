gnomAD <- read.csv("gnom.csv")
filter <- gnomAD[grep("^c.396", gnomAD$Consequence), ]

