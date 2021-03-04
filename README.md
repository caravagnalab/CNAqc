
# CNAqc <img src='man/figures/logo.png' align="right" height="139" />

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/caravagn/CNAqc.svg?branch=master)](https://travis-ci.org/caravagn/CNAqc)
[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![](https://img.shields.io/badge/Part%20of-evoverse-blue.svg)](https://caravagn.github.io/evoverse)
<!-- badges: end -->

`CNAqc` is a package that contains different methods to inspect the
quality, visualise and process absolute Copy Number Alteration (CNA)
calls generated from bulk sequencing of tumour samples.

Methods are available to integrate somatic mutation data of the same
tumour with CNA segments, tumour purity and ploidy assessment. The
package provides methods to estimate Cancer Cell Fractions (CCFs)
normalizing tumour mutation data for copy number states and tumour
purity. The package contains also statistical tests to identify patterns
of overfragmentation either in individual chromosome arms or at the
whole-genome level.

`CNAqc` is part of the `evoverse` set of [R
packages](https://caravagn.github.io/evoverse) to implement Cancer
Evolution analyses.

#### Reference

The paper that describes the CNAqc method is available as a preprint at
biorXiv:

*An automated quality checking tool for clonal copy number changes and
single nucleotide variant calls from high-resolution whole genome
sequencing* Jacob Househam, William CH Cross and Giulio Caravagna.
Preprint, 2021.

#### Help and support

## [![](https://img.shields.io/badge/GitHub%20Pages-https://caravagnalab.github.io/CNAqc/-yellow.svg)](https://caravagnalab.github.io/CNAqc)

### Installation

You can install the released version of `CNAqc` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("caravagnalab/CNAqc")
```

------------------------------------------------------------------------

#### Copyright and contacts

Giulio Caravagna. Cancer Data Science (CDS) Laboratory.

[![](https://img.shields.io/badge/Email-gcaravagn@gmail.com-steelblue.svg)](mailto:gcaravagn@gmail.com)
[![](https://img.shields.io/badge/CDS%20Lab%20Github-caravagnalab-seagreen.svg)](https://github.com/caravagnalab)
[![](https://img.shields.io/badge/CDS%20Lab%20webpage-https://www.caravagnalab.org/-red.svg)](https://www.caravagnalab.org/)
