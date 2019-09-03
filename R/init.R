#' Title
#'
#' @param snvs
#' @param cna
#' @param purity
#'
#' @return
#' @export
#'
#' @examples
init = function(snvs, cna, purity)
{
  input = prepare_input_data(snvs, cna, tumour_purity)

  fit = list()
  class(fit) <- "cnaqc"

  fit$snvs = input$snvs
  fit$cna = input$cna

  fit$n_snvs = nrow(snvs)
  fit$n_cna = nrow(cna)

  fit$n_cna_clonal = sum(fit$cna$CCF == 1)
  fit$n_cna_sbclonal = sum(fit$cna$CCF < 1)

  fit$n_karyotype = sort(table(fit$snvs$karyotype), decreasing = T)

  fit
}
