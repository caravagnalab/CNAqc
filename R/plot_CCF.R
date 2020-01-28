#' Plot the CCF estimates in the data.
#'
#' @description
#'
#' This function extracts the plots of the Cancer Cell Fractions (CCFs) estimates
#' from a `CNAqc` object, after they have been computed using
#' function `compute_CCF`. The plots - one per computed karyotype - are assembled
#' into a `ggpubr` figure.
#'
#' @param x An object of class \code{cnaqc}, where CCF have been computed using
#' function `compute_CCF`.
#'
#' @return A `ggpubr` figure of the CCF estiamtes.
#'
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc')
#' x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna, example_dataset_CNAqc$purity)
#'
#' x = compute_CCF(x)
#' plot_CCF(x)
plot_CCF = function(x)
{
  stopifnot(inherits(x, 'cnaqc'))

  with_CCF = all(!is.null(x$CCF_estimates))
  if(!with_CCF){
    stop("Input does not have CCF estimates, see ?compute_CCF to determine CCF values.")
  }

  ccf = lapply(names(x$CCF_estimates), plot_mutation_multiplicity_entropy, x = x)


  ggpubr::ggarrange(
    plotlist = ccf,
    ncol = 1,
    nrow = length(ccf)
  )
}
