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
  input = CNAqc:::prepare_input_data(snvs, cna, purity)

  fit = list()
  class(fit) <- "cnaqc"

  # Remove CNA segments with NA Major/minor
  na_allele_Major = sapply(input$cna$Major, is.na)
  na_allele_minor = sapply(input$cna$minor, is.na)
  na_allele = na_allele_Major | na_allele_minor
  if(length(na_allele) > 0)
  {
    cat(crayon::red("\n[CNAqc] The following segments have Major/ minor alleles in non-numeric format or NAs\n"))
    pio::pioDisp(input$cna[na_allele, , drop = FALSE])
    cat(crayon::red("These CNA segments will be removed.\n"))

    input$cna = input$cna[!na_allele, , drop = FALSE]
  }

  # TODO -- check NA chrom/from and to, and the same for SNVs

  fit$snvs = input$snvs
  fit$cna = input$cna

  # Counts data
  fit$n_snvs = nrow(snvs)
  fit$n_cna = nrow(cna)

  fit$n_cna_clonal = sum(fit$cna$CCF == 1)
  fit$n_cna_sbclonal = sum(fit$cna$CCF < 1)

  fit$n_karyotype = sort(table(fit$snvs$karyotype), decreasing = T)
  fit$purity = purity

  fit
}
