#' Print for class \code{'cnaqc'}.
#'
#' @param x An obj of class \code{'cnaqc'}.
#' @param ... Default S3 method parameter.
#'
#' @return Nothing.
#'
#' @export
#'
#' @import pio
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna,example_dataset_CNAqc$purity)
#'
#' print(x)
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

  pio::pioStr("\n     Purity ", paste0(x$purity  *100, ' cellularity'))
  pio::pioStr("\n Karyotypes ", paste0(x$n_karyotype, ' (', names(x$n_karyotype), ')', collapse = '; '), '\n')

  with_peaks = all(!is.null(x$peaks_analysis))

  if(with_peaks)
  {
    pio::pioStr("\n   Peaks QC ", with_peaks, ' ~ s = ', x$peaks_analysis$score, suffix = '\n')
    pio::pioDisp(x$peaks_analysis$matches)
  }
  else pio::pioStr("\n   Peaks QC ", FALSE, suffix = '\n')

}
