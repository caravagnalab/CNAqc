% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maf.R
\name{as_maftools_obj}
\alias{as_maftools_obj}
\title{Convert a CNAqc object to a maftools object.}
\usage{
as_maftools_obj(
  x,
  only_drivers = TRUE,
  CNA_genes = NULL,
  cross_reference = TRUE,
  assembly = TRUE
)
}
\arguments{
\item{x}{A CNAqc object with MAF annotations.}

\item{only_drivers}{If `TRUE`, only driver mutations are used, otherwised all.
When `TRUE`, if drivers are not annotated, an error is thrown.}

\item{CNA_genes}{The list of genes for ....}

\item{assembly}{If `TRUE`...}
}
\value{
a \code{"maftools"} object.
}
\description{
A single CNAqc object for which MAF annotations have been added using function
\code{\link{augment_with_maf}}, can be converted as an object for the Bioconductor
package [maftools](https://bioconductor.org/packages/release/bioc/html/maftools.html).

The export will use only driver mutations data or all the annotated mutations,
depending on the input flag `only_drivers`. Moreover, CNAqc copy number states
of a list of desired genes can also be augmented to the
[maftools](https://bioconductor.org/packages/release/bioc/html/maftools.html)
object.
}
\examples{
if(FALSE)
{
   # Create your CNAqc object (omissis here) from an original "file.vcf"
   x = init(mutations = ..., cna = ..., purity = ...)

   # Offline, create your MAF annotations as file "file_vcf.maf" from "file.vcf"
   # vcf2maf file.vcf .... file_vcf.maf

   # Import into R/CNAqc
   x = augment_with_maf(x, maf = "file_vcf.maf")

   # Extraction
   x \%>\% as_maftools_obj
}
}
\seealso{
function \code{\link{augment_with_maf}} to add MAF annotations to a
CNAqc object, to be used before running `as_maftools_obj`.
}
