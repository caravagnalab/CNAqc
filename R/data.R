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
#' @examples
#' data(example_PCAWG)
#' example_PCAWG
#'
"example_PCAWG"

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
#' @examples
#' data(fpr_test)
#' fpr_test
#'
"fpr_test"

#' List of Intogen driver genes per tumour type.
#'
#' @docType data
#'
#' @usage data(intogen_drivers)
#'
#' @format A tibble.
#'
#' @keywords datasets
#'
#' @references Francisco Martínez-Jiménez, Ferran Muiños, Inés Sentís, Jordi
#' Deu-Pons, Iker Reyes-Salazar, Claudia Arnedo-Pac, Loris Mularoni, Oriol Pich,
#' Jose Bonet, Hanna Kranas, Abel Gonzalez-Perez & Nuria Lopez-Bigas.
#' A compendium of mutational cancer driver genes Nat Rev Cancer (2020).;
#' doi:10.1038/s41568-020-0290-x
#' @examples
#' data(intogen_drivers)
#' intogen_drivers
#'
"intogen_drivers"

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

#' Coordinates for hg19 genes
#'
#' @docType data
#'
#' @usage data(gene_coordinates_hg19)
#'
#' @format A tibble that represents the coordinates for genes in the hg19 genome assembly,
#' reporting the gene symbol, from and to (chromosome range).
#'
#' @keywords datasets
#'
#' @references hg19
#'
#' @examples
#' data(gene_coordinates_hg19)
#' gene_coordinates_hg19
#'
"gene_coordinates_hg19"
