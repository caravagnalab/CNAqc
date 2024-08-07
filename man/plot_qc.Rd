% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_qc.R
\name{plot_qc}
\alias{plot_qc}
\title{Plot a summary of QC results.}
\usage{
plot_qc(x)
}
\arguments{
\item{x}{A CNAqc object.}
}
\value{
A \code{ggplot2} plot
}
\description{
Results from \code{analyze_peaks} and \code{compute_CCF} can be
visualised with this function. Compared to individual karyotypes fits available
with function \code{plot_peaks_analysis}, for instance, this function reports
summary pass/fail statistics for each analysis.
}
\examples{
data('example_dataset_CNAqc', package = 'CNAqc')
x = init(mutations = example_dataset_CNAqc$mutations, cna = example_dataset_CNAqc$cna, purity = example_dataset_CNAqc$purity)

x = analyze_peaks(x)
x = compute_CCF(x)

plot_qc(x)
}
