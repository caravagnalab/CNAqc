
# CNAqc <a href="caravagnalab.github.io/CNAqc"><img src="man/figures/logo.png" align="right" height="139" /></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/caravagnalab/CNAqc/workflows/R-CMD-check/badge.svg)](https://github.com/caravagnalab/CNAqc/actions)
[![pkgdown](https://github.com/caravagnalab/CNAqc/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/caravagnalab/CNAqc/actions/workflows/pkgdown.yaml)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://www.tidyverse.org/lifecycle/#stable)
<!-- badges: end -->

CNAqc is a package to quality control (QC) bulk cancer sequencing data.
Methods are available to , visualise and manipulate i) somatic mutation
data of both single-nucleotide variants and insertion-deletions, ii)
allele-specific Copy Number Alterations (CNAs) and iii) tumour purity
estimates. QC procedures in CNAqc can be used to validate copy number
segmentations against variant allele frequencies of somatic mutations;
QC scores can be used to rank alternative tumour segmentations and
purity/ ploidy estimates. The tool also provides an automatic copy
number calling pipeline that uses the Sequenza caller algorithm. CNAqc
provides also algorithms to phase mutation multiplicities against CNAs
and estimate Cancer Cell Fractions (CCFs) with their uncertainty. The
package contains also statistical tests to identify patterns of
over-fragmentation of chromosome arms (excessively short and numerous
DNA fragments) and perform various manipulation tasks for somatic tumour
data.

#### Citation

[![](https://img.shields.io/badge/doi-10.1101/2021.02.13.429885-red.svg)](https://doi.org/10.1101/2021.02.13.429885)

If you use `CNAqc`, please cite:

-   *Integrated quality control of allele-specific copy numbers,
    mutations and tumour purity from cancer whole genome sequencing
    assays* Jacob Househam, Riccardo Bergamin, Salvatore Milite, Nicola
    Calonaci, Alice Antonello, Marc J Williams, William CH Cross, Giulio
    Caravagna. [biorXiv
    2022](https://www.biorxiv.org/content/10.1101/2021.02.13.429885v3).

#### Help and support

## [![](https://img.shields.io/badge/GitHub%20Pages-https://caravagnalab.github.io/CNAqc/-yellow.svg)](https://caravagnalab.github.io/CNAqc)

### Installation

``` r
# install.packages("devtools")
devtools::install_github("caravagnalab/CNAqc")
```

------------------------------------------------------------------------

#### Copyright and contacts

Cancer Data Science (CDS) Laboratory, University of Trieste, Italy.

[![](https://img.shields.io/badge/CDS%20Lab%20Github-caravagnalab-seagreen.svg)](https://github.com/caravagnalab)
[![](https://img.shields.io/badge/CDS%20Lab%20webpage-https://www.caravagnalab.org/-red.svg)](https://www.caravagnalab.org/)
