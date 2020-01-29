#' Compute CCF values for the available mutations.
#'
#' @description
#'
#' This function provides the implementation of a set of entropy-related routines that can estimate
#' Cancer Cell Fraction values (CCFs) for the available mutations. The implemented routine is described
#' in the package vignette `"Computation of Cancer Cell Fractions"` that is available at the URL
#' \url{https://caravagn.github.io/CNAqc/articles/ccf_computation.html}. This function creates a field
#' `CCF_estimates` inside the returned object which contains both the estimated CCF values and the
#' plot of the report of this analysis.
#'
#' @param x An object of class \code{cnaqc}, created by the \code{init} function.
#' @param karyotypes The karyotypes to use, this package supports only \code{c('2:1', '2:0', '2:2')}.
#' @param entropy_quantile The entropy quantile used to determine the interval in which the CCF estimates
#' are not reliable because it is not possible to determine a precise value for the multiplicity of a
#' mutation from its allelic frequencies. See the package vignette `"Computation of Cancer Cell Fractions"`
#' that is available at the URL \url{https://caravagn.github.io/CNAqc/articles/ccf_computation.html} to
#' see the detailed meaning of this parameter.
#'
#' @seealso Getters function \code{CCF} and \code{plot_CCF}.
#' @return An object of class \code{cnaqc}, with CCF values available for extraction and plotting.
#' @export
#'
#' @examples
#'
#' data('example_dataset_CNAqc')
#' x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna, example_dataset_CNAqc$purity)
#'
#' x = compute_CCF(x, karyotypes = c('1:0', '1:1', '2:1', '2:0', '2:2'))
#' print(x)
#'
#' CCF(x)
#' plot_CCF(x)
compute_CCF = function(x, karyotypes = c('1:0', '1:1', '2:1', '2:0', '2:2'), entropy_quantile = .9)
{
  stopifnot(inherits(x, 'cnaqc'))

  nkaryotypes = x$n_karyotype[x$n_karyotype > 25]

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
        return(mutmult_single_copy(x, k))

      return(mutmult_two_copies(x, k, entropy_quantile = .9))
    })
  names(x$CCF_estimates) = karyotypes

  # Report some stats
  mutations = lapply(x$CCF_estimates , function(x) x$mutations)
  mutations = Reduce(bind_rows, mutations)

  cat('\n')
  cli::cli_h2("Summary CCF assignments. NA: not assignable with q = {.field {entropy_quantile}}")
  pioDisp(
    mutations %>%
      group_by(karyotype, mutation_multiplicity) %>%
      summarise(assignments = n()) %>%
      ungroup()
  )
  cat("Note: NA ~ mutations not confidently assignable with q =", entropy_quantile, '\n')

  x
}



