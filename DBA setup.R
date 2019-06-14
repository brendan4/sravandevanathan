
library("devtools")
devtools::install_github("klutometis/roxygen")
library(roxygen2)
create("DBA")
install("DBA")
setwd("./DBA")
document()
lsf.str("package:DBA")
