#' Coordinates for GRCh38 genes
#'
#' @docType data
#'
#' @usage data(gene_coordinates_GRCh38)
#'
#' @format A tibble that represents the coordinates for genes in the GRCh38 genome assembly,
#' reporting the gene symbol, from and to (chromosome range).
#'
#' @keywords datasets
#'
#' @references GRCh38
#'
#' @examples
#' data(gene_coordinates_GRCh38)
#' gene_coordinates_GRCh38
#'
"gene_coordinates_GRCh38"

#' Coordinates for hg19 chromosomes.
#'
#' @docType data
#'
#' @usage data(chr_coordinates_hg19)
#'
#' @format A tibble that represents the coordinates for hg19 genome assembly, reporting the chromosome label,
#' from and to (chromosome range), the length of the chromosome, the position (start and end) of the
#' centromers.
#'
#' @keywords datasets
#'
#' @references hg19
#'
#' @examples
#' data(chr_coordinates_hg19)
#' chr_coordinates_hg19
#'
"chr_coordinates_hg19"

#' Coordinates for GRCh38 chromosomes.
#'
#' @docType data
#'
#' @usage data(chr_coordinates_GRCh38)
#'
#' @format A tibble that represents the coordinates for GRCh38 genome assembly, reporting the chromosome label,
#' from and to (chromosome range), the length of the chromosome, the position (start and end) of the
#' centromers.
#'
#' @keywords datasets
#'
#' @references GRCh38
#'
#' @examples
#' data(chr_coordinates_GRCh38)
#' chr_coordinates_GRCh38
#'
"chr_coordinates_GRCh38"

#' Example CNAqc dataset.
#'
#' @docType data
#'
#' @usage data(example_dataset_CNAqc)
#'
#' @format A list of SNVs, allele-specific copy number alterations (CNAs) and purity value that
#' can be used with CNAqc. This tumour does not contain subclonal CNAs.
#'
#' @keywords datasets
#'
#' @examples
#' data(example_dataset_CNAqc)
#' example_dataset_CNAqc
#'
"example_dataset_CNAqc"

#' Example PCAWG tumour
#'
#' @docType data
#'
#' @usage data(example_PCAWG)
#'
#' @format A CNAqc object created for a PCAWG sample, with subclonal CNAs
#' called by Battenberg.
#'
#' @keywords datasets
#'
#' @references The ICGC/TCGA Pan-Cancer Analysis of Whole Genomes Consortium.
#' Pan-cancer analysis of whole genomes. Nature 578, 82–93 (2020).
#'  https://doi.org/10.1038/s41586-020-1969-6
#'
#' @examples
#' data(example_PCAWG)
#' example_PCAWG
#'
"example_PCAWG"

#' Driver mutations released with the DriverDBv3 database
#' \url{http://driverdb.tms.cmu.edu.tw/}
#'
#' @docType data
#'
#' @usage data(DriverDBv3)
#'
#' @format A tibble.
#'
#' @keywords datasets
#'
#' @references ??
#'
#' @examples
#' data(DriverDBv3)
#' DriverDBv3
#'
"DriverDBv3"

#' Data (simulation performance) from the trainig set to auto-tune epsilon.
#'
#' @docType data
#'
#' @usage data(fpr_test)
#'
#' @format A tibble.
#'
#' @keywords datasets
#'
#' @references
#'
#' @examples
#' data(fpr_test)
#' fpr_test
#'
"fpr_test"

#' Driver mutations released with the boost DM method the DriverDBv3 database
#' \url{https://www.intogen.org/boostdm/}
#'
#' @docType data
#'
#' @usage data(boostDM_drivers)
#'
#' @format A tibble.
#'
#' @keywords datasets
#'
#' @references Muiños, F., Martínez-Jiménez, F., Pich, O. et al.
#' In silico saturation mutagenesis of cancer genes. Nature 596, 428–432 (2021).
#' https://doi.org/10.1038/s41586-021-03771-1
#'
#' @examples
#' data(boostDM_drivers)
#' boostDM_drivers
#'
"boostDM_drivers"


