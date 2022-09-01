#' Compute CCF values for the available mutations.
#'
#' @description
#'
#' This function provides the implementation of a set of entropy-related routines that can estimate
#' Cancer Cell Fraction values (CCFs) for the available mutations. The implemented routine is described
#' in the package vignette `"Computation of Cancer Cell Fractions"` that is available at the URL
#' \url{https://caravagnalab.github.io/CNAqc/articles/ccf_computation.html}. This function creates a field
#' `CCF_estimates` inside the returned object which contains both the estimated CCF values and the
#' plot of the report of this analysis.
#'
#' @param x An object of class \code{cnaqc}, created by the \code{init} function.
#' @param karyotypes The karyotypes to use, this package supports only \code{c('2:1', '2:0', '2:2')}.
#' @param cutoff_QC_PASS Percentage of mutations that can be not-assigned (\code{NA}) in a karyotype.
#' If the karyotype has more than \code{cutoff_QC_PASS} percentage of non-assigned mutations then
#' the overall set of CCF calls is failed for the karyotype.
#' @param muts_per_karyotype Minimum number of mutations that are required to be mapped to a karyotype
#' in order to compute CCF values (default 25).
#' @param method Either \code{"ENTROPY"} or \code{"ROUGH"}, to reflect the two different algorithms
#' to compute CCF. See the package vignette to understand the differences across methods.
#'
#' @seealso Getters function \code{CCF} and \code{plot_CCF}.
#' @return An object of class \code{cnaqc}, with CCF values available for extraction and plotting.
#' @export
#'
#' @examples
#'
#' data('example_dataset_CNAqc')
#' x = init(example_dataset_CNAqc$mutations, example_dataset_CNAqc$cna, example_dataset_CNAqc$purity)
#'
#' x = compute_CCF(x, karyotypes = c('1:0', '1:1', '2:1', '2:0', '2:2'))
#' print(x)
#'
#' # Extract the values with these other functions
#' CCF(x)
#' plot_CCF(x)
compute_CCF = function(x,
                       karyotypes = c('1:0', '1:1', '2:0', '2:1', '2:2'),
                       muts_per_karyotype = 25,
                       cutoff_QC_PASS = 0.1,
                       method = 'ENTROPY')
{
  stopifnot(inherits(x, 'cnaqc'))
  stopifnot(method %in% c('ENTROPY', "ROUGH"))

  if(any(x$n_karyotype <= muts_per_karyotype)) warning("Some karyotypes have fewer than", muts_per_karyotype, 'and will not be analysed.')
  nkaryotypes = x$n_karyotype[x$n_karyotype > muts_per_karyotype]

  karyotypes = intersect(karyotypes, nkaryotypes %>% names)
  stopifnot(
    karyotypes %in% c('1:0', '1:1', '2:1', '2:0', '2:2')
  )

  # Compute mutation multiplicity
  x$CCF_estimates = lapply(
    karyotypes,
    function(k)
    {
      if(k %in% c('1:0', '1:1'))
        return(suppressWarnings(CNAqc:::mutmult_single_copy(x, k)))

      if(method == "ENTROPY")
        return(suppressWarnings(CNAqc:::mutmult_two_copies_entropy(x, k)))
      else
        return(suppressWarnings(CNAqc:::mutmult_two_copies_rough(x, k)))
    })
  names(x$CCF_estimates) = karyotypes

  # Check if there is any null (errors), and remove it
  null_entries = sapply(x$CCF_estimates, function(x) all(is.null(x)))
  x$CCF_estimates = x$CCF_estimates[!null_entries]

  # On extreme cases where there is NO CCF available, we just return x
  if(length(x$CCF_estimates) == 0) {
    x$CCF_estimates = NULL
    return(x)
  }

  # Report some stats
  mutations = lapply(x$CCF_estimates , function(x) x$mutations)
  mutations = Reduce(dplyr::bind_rows, mutations)

  # pioDisp(
  #   mutations %>%
  #     dplyr::group_by(karyotype, mutation_multiplicity) %>%
  #     dplyr::summarise(assignments = n()) %>%
  #     dplyr::ungroup()
  # )

  # QC the findings
  N = mutations %>%
    dplyr::group_by(karyotype) %>%
    dplyr::summarise(N = n()) %>%
    dplyr::ungroup()

  NA_N = mutations %>%
    dplyr::group_by(karyotype, mutation_multiplicity) %>%
    dplyr::summarise(Unknown = n()) %>%
    dplyr::filter(is.na(mutation_multiplicity)) %>%
    dplyr::select(-mutation_multiplicity) %>%
    dplyr::ungroup()

  QC_table = N %>%
    dplyr::full_join(NA_N, by = 'karyotype') %>%
    dplyr::mutate(
      Unknown = ifelse(is.na(Unknown), 0, Unknown),
      p_Unkown = Unknown/N,
      QC = ifelse(p_Unkown < cutoff_QC_PASS, "PASS", "FAIL"),
      method = method
    )

  if(any(QC_table$QC == "FAIL"))
  {
    cat('\n')
    cli::cli_h2("Summary CCF assignments. (>{.field {cutoff_QC_PASS*100}%} NAs: not assignable with confidence)")
    pioDisp(QC_table)
  }

  for(k in QC_table$karyotype)
    x$CCF_estimates[[k]]$QC_table = QC_table %>% dplyr::filter(karyotype == !!k)

  x
}



