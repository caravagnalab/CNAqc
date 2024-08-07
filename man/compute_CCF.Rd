% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_CCF.R
\name{compute_CCF}
\alias{compute_CCF}
\title{Compute CCF values.}
\usage{
compute_CCF(
  x,
  karyotypes = c("1:0", "1:1", "2:0", "2:1", "2:2"),
  muts_per_karyotype = 25,
  cutoff_QC_PASS = 0.1,
  method = "ENTROPY"
)
}
\arguments{
\item{x}{A CNAqc object.}

\item{karyotypes}{The karyotypes to use, this package supports only clonal simple CNAs.}

\item{muts_per_karyotype}{Minimum number of mutations that are required to be mapped to a karyotype
in order to compute CCF values (default 25).}

\item{cutoff_QC_PASS}{For the entropy-based method, percentage of mutations that
can be not-assigned (\code{NA}) in a karyotype. If the karyotype has more than
\code{cutoff_QC_PASS} percentage of non-assigned mutations, then the overall set of CCFs
is failed for the karyotype.}

\item{method}{Either \code{"ENTROPY"} (default) or \code{"ROUGH"}, to reflect the two different algorithms
to compute CCF.}
}
\value{
A CNAqc object with CCF values.
}
\description{
With this function, CNAqc can compute a CCF per mutation upon “phasing” the multiplicity
for every input mutation. Phasing is the task of computing the number of copies of a
mutation mapping in a certain copy number segment; this task is a difficult, and can lead to erroneous CCF estimates.

CNAqc computes CCFs for simple clonal CNA segments, offering two algorithms to phase
mutations directly from VAFs.

* Entropy based method. The entropy-based approach will flag mutations for which
 we cannot phase multiplicity by VAFs with certainty; the CCFs of these mutations
 should be manually controlled and, unless necessary, discarded. To this aid, a QC pass
 is assigned with less than a certain percentage of mutations have uncertain CCFs. The model uses the
 entropy of a VAF mixture with two Binomial distributions to detect mutations happened
 before, and after aneuploidy. Assigning multiplicities is difficult at the crossing
 of the two densities where mutations could have multiplicity 1 or 2. If mistaken,
 these mutations can determine aritficial peaks in the CCF distribution and compromise
 downstream subclonal deconvolution.

 * Hard-cut based method. A method is available to compute CCFs regardless of the
 entropy. From the 2-class Binomial mixture, CNAqc uses the means of the Binomial
 parameters to determine a hard split of the data. Since there are no NA assignments,
 the computation is always scored PASS for QC purposes; for this reason this computation is more “rough” than the one based on entropy.

Like for other analyses This function creates a field `CCF_estimates` inside
the returned object which contains the estimated CCFs.
}
\examples{
data('example_dataset_CNAqc')
x = init(mutations = example_dataset_CNAqc$mutations, cna =example_dataset_CNAqc$cna, purity = example_dataset_CNAqc$purity)

x = compute_CCF(x, karyotypes = c('1:0', '1:1', '2:1', '2:0', '2:2'))
print(x)

# Extract the values with these other functions
CCF(x)
plot_CCF(x)
}
\seealso{
Getters function \code{CCF} and \code{plot_CCF}.
}
