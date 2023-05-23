#' Extract CCF estimates.
#'
#' @description
#'
#' This function extracts Cancer Cell Fractions (CCFs) estimates
#' from a CNAqc object, after they have been computed using
#' function `compute_CCF`. The estimates are pooled
#' across the used karyotypes, and extracted from the
#' `CCF_estimates` field of the input object.
#'
#' @param x A CNAqc object.
#'
#' @return A mutations tibble with columns for mutation multiplicity and CCFs.
#'
#' @export
#'
#' @examples
#' data("example_PCAWG", package = 'CNAqc')
#' CCF(example_PCAWG)
CCF = function(x)
{
  stopifnot(inherits(x, 'cnaqc'))

  with_CCF = all(!is.null(x$CCF_estimates))
  if(!with_CCF){
    stop("Input does not have CCF estimates, see ?compute_CCF to determine CCF values.")
  }

  mutations = lapply(x$CCF_estimates , function(x) x$mutations)
  mutations = Reduce(bind_rows, mutations) %>% mutate(cna = "clonal")

  return(mutations)
}

#' Extract mutations.
#'
#' @description Getter to obtain mutation calls from an object.
#'
#' @param x A CNAqc object.
#' @param cna \code{"clonal"} for clonal CNAs, \code{"subclonal"} for subclonal CNAs.
#' @param type \code{"SNV"} for single-nucleotide variants, \code{"indel"} for insertion-deletions.
#'
#' @return A tibble with the data.
#' @export
#'
#' @examples
#' data("example_PCAWG", package = 'CNAqc')
#' Mutations(example_PCAWG)
Mutations = function(x, cna = c("clonal", "subclonal"), type = c("SNV", "indel"))
{
  stopifnot(inherits(x, 'cnaqc'))

  clonal = NULL
  if("clonal" %in% cna)
    clonal = x$mutations %>% mutate(cna = 'clonal')

  subclonal = NULL
  if(("subclonal" %in% cna) & x$has_subclonal_CNA)
    subclonal = x$cna_subclonal$mutations %>% Reduce(f = bind_rows) %>% mutate(cna = 'subclonal')

  mutations = bind_rows(clonal, subclonal) %>%
    dplyr::select(chr, from, to, ref, alt, NV, DP, VAF, everything()) %>%
    filter(type %in% !!type)

  if((mutations %>% nrow())== 0) cli::cli_alert_danger("No mutations with these parameters: CNA {.field {cna}}, type {.field {type}}.")

  return(mutations)
}

#' Extract CNAs.
#'
#' @description Getter to obtain copy number calls from an object.
#'
#' @param x A CNAqc object.
#' @param type \code{"clonal"} for clonal CNAs, \code{"subclonal"} for subclonal CNAs.
#'
#' @return A tibble with the data.
#' @export
#'
#' @examples
#' data("example_PCAWG", package = 'CNAqc')
#' CNA(example_PCAWG)
CNA = function(x, type = c("clonal", "subclonal"))
{
  stopifnot(inherits(x, 'cnaqc'))

  clonal = NULL
  if("clonal" %in% type)
    clonal = x$cna

  subclonal = NULL
  if(("subclonal" %in% type) & x$has_subclonal_CNA)
    subclonal = x$cna_subclonal

  cna = dplyr::bind_rows(clonal, subclonal) %>%
    dplyr::select(chr, from, to, starts_with('Major'), starts_with('minor'), CCF, everything())

  if((cna %>% nrow())== 0) cli::cli_alert_danger("No CNAs with these parameters: {.field {cna}}.")

  return(cna)
}

#' Returns percentage of passed segments
#'
#' @description Getter to the percentage of genome with pass value
#'
#' @param x A CNAqc object.
#'
#' @return A tibble with the data.
#' @export
#'
#' @examples
#' data("example_PCAWG", package = 'CNAqc')
#' get_PASS_percentage(example_PCAWG)
get_PASS_percentage <- function(x) {
  stopifnot(inherits(x, 'cnaqc'))
  
  tb_cna <- x$cna %>% dplyr::group_by(QC_PASS) %>% dplyr::summarize(Length = sum(length) / sum(x$cna$length), N_segments = n()/ nrow(x$cna))
  tb_muts <- x$mutations %>% dplyr::group_by(QC_PASS) %>% dplyr::summarize(N_mutations = n()/ nrow(x$mutations))
  
  return(dplyr::full_join(tb_cna, tb_muts))
}




