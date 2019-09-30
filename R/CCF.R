#' Return the available mutation CCF estimates in the data.
#'
#' @description
#'
#' This function extracts Cancer Cell Fractions (CCFs) estimates
#' from a `CNAqc` object, after they have been computed using
#' function `compute_CCF`. The estimates are pooled
#' across the used karyotypes, and extracted from the
#' `CCF_estimates` field of the input object.
#'
#' @param x An object of class \code{cnaqc}, where CCF have been computed using
#' function `compute_CCF`.
#'
#' @return A tibble with the mutations with the new columns for mutation multiplicity
#' and CCF values.
#'
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc')
#' x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna, example_dataset_CNAqc$purity)
#'
#' x = compute_CCF(x)
#' CCF(x)
CCF = function(x)
{
  stopifnot(inherits(x, 'cnaqc'))

  with_CCF = all(!is.null(x$CCF_estimates))
  if(!with_CCF){
    stop("Input does not have CCF estimates, see ?compute_CCF to determine CCF values.")
  }

  mutations = lapply(x$CCF_estimates , function(x) x$mutations)
  mutations = Reduce(bind_rows, mutations)

  return(mutations)
}
