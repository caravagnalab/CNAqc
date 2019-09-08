#' Print for an object of class \code{'cnaqc'}.
#'
#' @param x An obj of class \code{'dbpmm'}.
#' @param ...
#'
#' @return nothing.
#' @export
#' @import crayon
#'
#' @examples
print.cnaqc = function(x, ...)
{
  stopifnot(inherits(x, "cnaqc"))

  # pio::pioHdr("[CNAqc - CNA Quality Check]")

  # Mapping mutations
  pio::pioStr("      CNAqc ",
              'n =',
              x$n_snvs,
              "mutations for",
              x$n_cna,
              paste0("CNA segments (", x$n_cna_clonal, " clonal, ", x$n_cna_sbclonal, " subclonal)"))

  pio::pioStr("\n   Purity ", paste0(x$purity  *100, ' cellularity'))


  pio::pioStr("\n Karyotypes ", paste0(x$n_karyotype, ' (', names(x$n_karyotype), ')', collapse = '; '), '\n')
}
