% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_Battenberg.R
\name{parse_Battenberg}
\alias{parse_Battenberg}
\title{Parse Battenberg calls.}
\usage{
parse_Battenberg(x)
}
\arguments{
\item{x}{The input dataframe, required to have the following columns: \code{chr},
\code{from}, \code{to}, \code{battenberg_nMaj1_A}, \code{battenberg_nMin1_A},
\code{battenberg_nMaj2_A}, \code{battenberg_nMin2_A}, \code{battenberg_frac1_A},
and \code{battenberg_frac2_A}.}
}
\value{
Two tibbles in CNAqc-ready format.
}
\description{
Returns two tibbles in CNAqc-ready format, one with clonal calls
and the other with subclonal calls extracted by parsing the Battenberg format.
This function can be used to process raw date released by PCAWG.
}
\examples{
\dontrun{
# Load some CSV results from Battenberg
x = read.csv(....) \%>\% parse_Battenberg()

# Work with clonal calls (omitting mutations and other parameters here)
x = init(cna = x$clonal, ...)


# work with subclonal calls (omitting mutations and other parameters here)
x = analyze_peaks_subclonal(cna = x$subclonal, ...)
}
}
