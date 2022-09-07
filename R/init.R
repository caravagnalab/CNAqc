#' Creates a CNAqc object.
#'
#' @description Creates a CNAqc object from input
#' mutations, allele-specifc copy numbers and a tumour purity value.
#' The resulting object retains the input mutations that map on top
#' of the copy number segments, and allows for the computation
#' of the QC metrics available in the CNAqc package.
#'
#' Genomic coordinates in relative (per-chromosome) format are
#' transformed into absolute coordinates by means of a reference
#' genome providing the length of each chromosome. CNAqc supports
#' `hg19`/`GRCh37` and `hg38`/`GRCh38` references, which are embedded
#' into the package as `CNAqc::chr_coordinates_hg19` and
#' `CNAqc::chr_coordinates_GRCh38`. An abitrary reference can also be
#' provided it is stored in am equivalent format.
#'
#' @param mutations A dataframe of mutations with the following fields:
#' * `chr` chromosome name, e.g., \code{"chr3"}, \code{"chr8"}, \code{"chrX"}, ...
#' * `from` where the mutation start, an integer number
#' * `to` where the mutation ends, an integer number
#' * `ref` reference allele, e.g., \code{"A"}, \code{"ACC"}, \code{"AGA"}, ...
#' * `alt` alternative allele, e.g., \code{"A"}, \code{"ACC"}, \code{"AGA"}, ...
#' * `DP` sequencing depth at the locus, an integer number
#' * `NV` number of reads with the variant at the locus, an integer number
#' * `VAF` variant allele frequency (VAF), defined as `NV/DP`, at the locus, a real number in [0,1]
#'
#' @param snvs Deprecated parameter.
#'
#' @param cna A dataframe of allele-specific copy number with the following fields:
#' * `chr` chromosome name, e.g., \code{"chr3"}, \code{"chr8"}, \code{"chrX"}, ...
#' * `from` where the segment start, an integer number
#' * `to` where the segment ends, an integer number
#' * `Major` for the number of copies of the major allele (or A-allele), an integer number
#' * `minor` for the number of copies of the major allele (or B-allele), an integer number
#' * `CCF` an optional cancer cell fraction (CCF) column distinguishing clonal and subclonal segments, a real number in [0,1]
#' * `Major_2` optional for the number of copies of the major allele (or A-allele) in the second clone if present, an integer number
#' * `minor_2` optional for the number of copies of the major allele (or B-allele) in the second clone  if present, an integer number
#'
#' If the `CCF` value is present and equal to 1, a segment is considered clonal, otherwise
#' subclonal. If a segment is subclonal
#' * the columns `Major` and `minor` are interpreted as those for a subclone with proportion equal to the `CCF` value;
#' * the columns `Major_2` and `minor_2` are interpreted as those for a second subclone with proportion equal to the `1 - CCF` value;
#'
#' @param purity Value in between `0` and `1` to represent the proportion
#' of actual tumour content (sometimes called "cellularity").
#'
#' @param ref A key word for the used reference coordinate system. CNAqc supports
#' `hg19`/`GRCh37` and `hg38`/`GRCh38` references, which are embedded
#' into the package as `CNAqc::chr_coordinates_hg19` and
#' `CNAqc::chr_coordinates_GRCh38`. An abitrary reference can also be
#' provided if `ref` is a dataframe in the same format as `CNAqc::chr_coordinates_hg19`
#' or `CNAqc::chr_coordinates_GRCh38`. The default reference is `GRCh38`.
#'
#' @return A CNAqc object of class `cnaqc`, with S3 methods for printing,
#' plotting and analyzing data.
#'
#' @export
#'
#' @import crayon
#' @import vcfR
#' @import clisymbols
#' @import easypar
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' print(example_dataset_CNAqc)
#'
#' x = init(mutations = example_dataset_CNAqc$mutations, cna = example_dataset_CNAqc$cna, purity = example_dataset_CNAqc$purity)
#' print(x)
init = function(mutations, snvs = NULL, cna, purity, ref = "GRCh38")
{
  cli::cli_h1("CNAqc - CNA Quality Check")
  cat('\n')

  if(!is.null(snvs))
  {
    if(!is.null(mutations))
    {
      cli::boxx("Parameter `snvs` has been deprecated, cannot use it with `mutations`", col = 'red', margin = 3)

      cli::cli_abort("Avoid using `snvs` if you use `mutations")
    }
    else
    {
      cli::boxx("Parameter `snvs` has been deprecated, using it as `mutations",
                col = 'red',
                margin = 3) %>%
        cat('\n')

      mutations = snvs
    }

  }

  # Output
  fit = list()
  class(fit) <- "cnaqc"

  # Reference genome
  fit$reference_genome = ref
  cli::cli_alert_info("Using reference genome coordinates for: {.field {ref}}.")

  # Parse input
  input = prepare_input_data(mutations, cna, purity)

  # Remove CNA segments with NA Major/minor
  # na_allele_Major = sapply(input$cna_clonal$Major, is.na)
  # na_allele_minor = sapply(input$cna_clonal$minor, is.na)
  # na_allele = na_allele_Major | na_allele_minor
  #
  # discarded_cna = input$cna[na_allele, , drop = FALSE]
  # if(nrow(discarded_cna) > 0)
  # {
  #   cat(crayon::red("\n[CNAqc] CNA calls: the following segments have Major/ minor alleles in non-numeric format or NAs, and will be removed\n"))
  #   pio::pioDisp(discarded_cna)
  #
  #   input$cna = input$cna[!na_allele, , drop = FALSE]
  # }

  fit$mutations = input$mutations
  fit$cna = input$cna_clonal %>%
    dplyr::left_join(input$tab, by = 'segment_id') %>% as_tibble()
  fit$cna_subclonal = input$cna_subclonal %>% as_tibble()
  fit$has_subclonal_CNA = !all(is.null(input$cna_subclonal))

  # Counts data
  if(!fit$has_subclonal_CNA)
    fit$n_mutations = nrow(fit$mutations)
  else
    fit$n_mutations = nrow(fit$mutations)  + sapply(fit$cna_subclonal$mutations, nrow) %>% sum()

  fit$n_cna_clonal = nrow(fit$cna)
  fit$n_cna_subclonal = ifelse(is.null(fit$cna_subclonal), 0, nrow(fit$cna_subclonal))
  fit$n_cna = fit$n_cna_clonal + fit$n_cna_subclonal

  fit$n_karyotype = sort(table(fit$mutations$karyotype), decreasing = T)
  fit$purity = purity

  # Segments length (clonal)
  genome_segs_length = fit$cna %>%
    dplyr::group_by(Major, minor) %>%
    dplyr::summarise(L = sum(length), .groups = 'drop') %>%
    dplyr::mutate(karyotype = paste0(Major, ':', minor)) %>%
    dplyr::arrange(desc(L))

  fit$l_karyotype = genome_segs_length$L
  names(fit$l_karyotype ) = genome_segs_length$karyotype

  tab_ploidy = fit$cna %>%
    dplyr::group_by(minor, Major) %>%
    dplyr::summarise(n = sum(length)) %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::mutate(karyotype = paste0(Major, ':', minor)) %>%
    dplyr::ungroup()

  fit$ploidy = as.numeric(tab_ploidy$minor[1]) + as.numeric(tab_ploidy$Major[1])
  fit$most_prevalent_karyotype = paste0(tab_ploidy$Major[1], ':',  tab_ploidy$minor[1])
  fit$basepairs_by_karyotype = tab_ploidy

  fit$most_mutations_karyotype = names(fit$n_karyotype)[1]

  fit
}

check_custom_reference = function(x)
{
  required = CNAqc::chr_coordinates_hg19 %>% colnames()

  if(!all(colnames(x) %in% required))
  {
    cli::boxx("Problems with your custom reference", background_col = 'red', col = 'white') %>% cat('\n')

    cli::cli_alert_danger("This type of dataframe should be used, but you are missing columns")
    CNAqc::chr_coordinates_hg19 %>% print()

    cli::cli_abort("Cannot use your custom refernece")
  }
}

