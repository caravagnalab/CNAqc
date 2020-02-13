#' Create a CNAqc object.
#'
#' @description Creates a CNAqc object from a set of
#' mutation and CNA calls, and a tumour purity value.
#' The resulting object contains the mutations mapped
#' into CNA segments, and allows for the computation
#' of the QC metrics available in the CNAqc package.
#' CNAqc supports `hg19` coordinates.
#'
#' @param snvs A data.frame of mutations with the following fields
#' available: `chr`, `from`, `to` as `hg19` chromosome coordinates.
#' `ref` and `alt` for the reference and alternative alleles, and
#' `DP`, `NV` and `VAF` for the depth (total number of reads with both
#' reference and alternative), the number of variants with the alternative
#' allele and the allele frequency (VAF = NV/DP). Chromosome names must be
#' in the format `chr1`, `chr2`, etc.; alleles should be characters and
#' all other fields numeric.
#' @param cna A data.frame of CNA with the following fields
#' available: `chr`, `from`, `to` as `hg19` chromosome coordinates.
#' `Major` and `minor` for the number of copies of the major allele,
#' and the minor (B-allele). A `CCF` column can be used as a real-value
#' in between 0 and 1 to represent Cancer Cell Fractions; if this is not
#' available `CCF = 1` is set and all calls will refer to clonal segments.
#' Otherwise, segments with `CCF<1` would be considered subclonal CNAs.
#' @param purity Value in between `0` and `1` to represent the proportion
#' of actual tumour content (sometimes called "cellularity").
#' @param ref The reference genome (either "hg19" or "GRCh38"); the default is "GRCh38".
#'
#' @return A CNAqc object of class `cnaqc`, with S3 methods for printing,
#' plotting and analyzing data.
#'
#' @export
#'
#' @import tidyverse
#' @import pio
#' @import crayon
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' print(example_dataset_CNAqc)
#'
#' x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna, example_dataset_CNAqc$purity)
#' print(x)
init = function(snvs, cna, purity, ref = "GRCh38")
{
  pio::pioHdr("CNAqc - CNA Quality Check")
  cat('\n')
  
  # Output
  fit = list()
  class(fit) <- "cnaqc"
  
  # Reference genome
  fit$reference_genome = ref
  cli::cli_alert_info("Using reference genome coordinates for: {.field {ref}}.")
  
  # Partse input
  input = CNAqc:::prepare_input_data(snvs, cna, purity)
  
  # Remove CNA segments with NA Major/minor
  na_allele_Major = sapply(input$cna$Major, is.na)
  na_allele_minor = sapply(input$cna$minor, is.na)
  na_allele = na_allele_Major | na_allele_minor

  discarded_cna = input$cna[na_allele, , drop = FALSE]

  if(nrow(discarded_cna) > 0)
  {
    cat(crayon::red("\n[CNAqc] CNA calls: the following segments have Major/ minor alleles in non-numeric format or NAs, and will be removed\n"))
    pio::pioDisp(discarded_cna)

    input$cna = input$cna[!na_allele, , drop = FALSE]
  }

  # TODO -- check NA chrom/from and to, and the same for SNVs

  fit$snvs = input$snvs
  fit$cna = input$cna %>%
    left_join(input$tab, by = 'segment_id')

  # Counts data
  fit$n_snvs = nrow(fit$snvs)
  fit$n_cna = nrow(fit$cna)

  fit$n_cna_clonal = sum(fit$cna$CCF == 1)
  fit$n_cna_sbclonal = sum(fit$cna$CCF < 1)

  fit$n_karyotype = sort(table(fit$snvs$karyotype), decreasing = T)
  fit$purity = purity


  tab_ploidy = fit$cna %>%
    group_by(minor, Major) %>%
    summarise(n = sum(length)) %>%
    arrange(desc(n))

  fit$ploidy = tab_ploidy$minor[1] + tab_ploidy$Major[1]

  fit
}
