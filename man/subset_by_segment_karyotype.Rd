% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subsample.R
\name{subset_by_segment_karyotype}
\alias{subset_by_segment_karyotype}
\title{Subset calls by segments' karyotype.}
\usage{
subset_by_segment_karyotype(x, karyotypes)
}
\arguments{
\item{x}{An object of class \code{cnaqc}, created by the \code{init} function.}

\item{karyotype}{A list of karyotype ids in \code{"Major:minor"} notation
(e.g., \code{"1:1", "2,1", ...}) that will be retained.}
}
\value{
An object of class \code{cnaqc}, created by the \code{init} function.
}
\description{
For an object already created, it subsets the calls to those that involve
segments with a certain karyotype.
}
\examples{
data('example_dataset_CNAqc', package = 'CNAqc')
x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna,example_dataset_CNAqc$purity)

subset_by_segment_karyotype(x, '2:2')
subset_by_segment_karyotype(x, '2:1')
}