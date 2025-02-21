% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean.r
\name{clean.environment}
\alias{clean.environment}
\title{Clean environment}
\usage{
clean.environment(keep = NULL)
}
\arguments{
\item{keep}{all other objects in environment to keep. Pass objects in as a string (pad with quotes). 
Pass multiple objects in as a vector. If object in keep not found in environment a warning will be issued
and will terminate the function before removing anything}
}
\value{
Returns nothing but a clean environment
}
\description{
Cleans global enviroment after a long day of working. Always keeps expressed.trans and expressed.genes
}
\examples{
clean.environment(keep = c("pheno", "expressed.genesN", "gene.results", "bg_filt", "PCA.filt", "filt.set", "filt.names"))
clean.environment()
clean.enivroment(keep = "bg_filt")
}
