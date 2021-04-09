% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drivers.R
\name{get_drivers}
\alias{get_drivers}
\title{Extract driver data.}
\usage{
get_drivers(x, chromosomes = paste0("chr", c(1:22, "X", "Y")), which = "VAF")
}
\arguments{
\item{x}{An object of class \code{cnaqc}, created by the \code{init} function
that can have drivers annotated.}

\item{chromosomes}{Which chromosome to subset}

\item{which}{A keyword for \code{"VAF"} or \code{"CCF"}.}
}
\value{
A tibble
}
\description{
Returns a tibble with the drivers annotated on certain chromsomes,
extracting VAF or CCF values (if computed). If drivers are not
annotated \code{NULL} is returned.
}
\examples{
data('example_dataset_CNAqc', package = 'CNAqc')
x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna,example_dataset_CNAqc$purity)

get_drivers(x)
}