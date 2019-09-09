# CNAqc

CNAqc - Copy Number Alteration (CNA) quality check - is a package
to provide a set of metrics to assess the quality of CNA calls.

The package exploits peak detection methods to determine 

The package provides statistical measures to quantify the concordance 
between mutation and Copy Number calls, exploitinng the relation 
between allelic imbalance in CNA segements and allelic frequencies 
of somatic mutations. Quantitative metrics and plots for data 
exploration are available.

Statistical Model
-----------------
CNAqc implements a linear score where the relative size of each of a
set of karyotypes is used to weight the offset between an estimated
peak in the data, and its expectation. The expectations are determined
by standard CNA computations accouting for normal plodiy, tumour purity
and tumor ploidy. The peaks are determined after a KDE of the data, run
through a dedicated peak-detection package.

Motivation
----------
With this package it is easy to visually assess the concordance of
somatic mutation calls, and the CNA segements where they map. Quantitative
measures can be used to suggest adjustemnts of the current estimates
as modifications of the input purity.

***
**Author:** [Giulio Caravagna](https://sites.google.com/site/giuliocaravagna/), _Institute of Cancer Research, UK_.

**Contact:** [[@gcaravagna](https://twitter.com/gcaravagna); [giulio.caravagna@icr.ac.uk](mailto:giulio.caravagna@icr.ac.uk)]






