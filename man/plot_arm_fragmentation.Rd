% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_arm_fragmentation.R
\name{plot_arm_fragmentation}
\alias{plot_arm_fragmentation}
\title{Plot the arm level fragmentation test.}
\usage{
plot_arm_fragmentation(x, zoom = 0)
}
\arguments{
\item{x}{A CNAqc object.}

\item{zoom}{Number of maximum zoom panels to show in the bottom of
the figure. By default 0.}
}
\value{
A \code{ggpubr} figure
}
\description{
This function produces a arm-level report for the fragmentation test, with:

* a scatter of the counts per arm, with scaled the p-values;
* a jump statistics per arm, $J$,  the sum of the variation in total copy number
 profiles, evaluated among each pair of contiguous segments.

 Significantly overfragmented arms with high $J$ have a “scattered” copy number
  profile. Those with low $J$ are more uniform, as they show little no copy
  number change, and can be possibly smoothed (see below).
}
\examples{
data('example_dataset_CNAqc', package = 'CNAqc')
x = init(mutations = example_dataset_CNAqc$mutations, cna = example_dataset_CNAqc$cna, purity = example_dataset_CNAqc$purity)

x = detect_arm_overfragmentation(x, genome_percentage_cutoff = .2)
plot_arm_fragmentation(x)
}
