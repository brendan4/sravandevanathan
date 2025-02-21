% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prettygenes.R
\name{pretty.gene.name}
\alias{pretty.gene.name}
\title{Pretty gene names}
\usage{
\method{pretty}{gene.name}(data.set, as.row.names = FALSE,
  remove.dups = FALSE, col = NULL)
}
\arguments{
\item{data.set}{a data.frame with genes to be made pretty}

\item{remove.dups}{logical statment indicating duplicated rows in pretty col should be removed. 
Will always tigger a warning as information across the row is lost. Primarly used to set row names if dups found.}

\item{col}{an argument only used if the rownames are not where gene names are stored. Should be an integer representing
which col should be used.}

\item{as.rom.names}{logical value indicating if the pretty col should turned into row names. 
Warning will be triggered if duplicates found and will not execute rownames switch.}
}
\value{
returns a data.frame with a pretty col containing gene names or a data.frame with pretty rownames. 
Rownames will only be added if no duplicates found or remove duplicates specified.
}
\description{
Takes a data.frame with genes append to some unique identifer. EX: RPL11_21321_342342
returns RPL11
}
