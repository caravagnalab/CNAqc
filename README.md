
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CNAqc

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/caravagn/CNAqc.svg?branch=master)](https://travis-ci.org/caravagn/CNAqc)
<!-- badges: end -->

CNAqc is a package to provide a set of metrics to assess the quality of 
Copy Number Alteration (CNA) calls.

The package provides statistical measures to quantify the concordance 
between mutation and Copy Number calls, exploiting allelic imbalance in absolute CNA segments and allelic frequencies 
of somatic mutations. Quantitative metrics and plots for data 
exploration and quality check are available, allowing for a intuitive
assessment of the quality of calls. Quantitative measures can also be 
used to suggest adjustemnts of the current purity estimates to increase
the quality of CNA calls.

#### Statistical model

CNAqc implements a linear score where the relative size of each of a
set of karyotypes is used to weight the offset between an estimated
peak in the data, and its expectation. The expectations are determined
by standard CNA computations accouting for normal plodiy, tumour purity
and tumor ploidy. The peaks are determined after a KDE of the data, run
through a dedicated peak-detection package.

-----

### Installation

You can install the released version of CNAqc from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("caravagn/CNAqc")
```

-----

#### Copyright and contacts

Giulio Caravagna, PhD. _Institute of Cancer Research, London, UK_.

* Personal webpage: [https://bit.ly/2kc9E6Y](https://sites.google.com/site/giuliocaravagna/), 
* Email address: [giulio.caravagna@icr.ac.uk](mailto:giulio.caravagna@icr.ac.uk) and [gcaravagn@gmail.com](mailto:gcaravagn@gmail.com)
* Twitter feed: [@gcaravagna](https://twitter.com/gcaravagna)
* Other works: [caravagn](https://github.com/caravagn)
