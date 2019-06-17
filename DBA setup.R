
#intial set up
library("devtools")
devtools::install_github("klutometis/roxygen")
library(roxygen2)
create("DBA")
install("DBA")

# making docments 
setwd("./DBA")
document()


# listing properties pf package 
lsf.str("package:DBA")

# saving data in the package
devtools::use_data(expressed.trans, DBA)
devtools::use_data(expressed.genes, DBA)

# loading data
data("expressed.genes", package = "DBA") 

# undating documentation 
devtools::document()
?filter.genes
?var.samples
?mergeTables

