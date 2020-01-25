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

  cli::cli_rule(
    paste(
      crayon::bgYellow(crayon::black("[ CNAqc ] ")),
      'n = {.value {x$n_snvs}} mutations in {.field {x$n_cna}} segments ({.value {x$n_cna_clonal}} clonal + {.value {x$n_cna_sbclonal}} subclonal)'
    )
  )

  # cli::cli_alert_info(paste0(" CNA segments: ", x$n_cna_clonal, " clonal, ", x$n_cna_sbclonal, " subclonal."))
  cli::cli_alert_info(paste0("Purity: ", paste0(x$purity  *100, '% ~ Ploidy: ', x$ploidy, '.')))

  cli::cli_alert_info(paste0("Mutation mapping (head): ", paste0(head(x$n_karyotype), ' (',
                                                        names(head(x$n_karyotype)), ')', collapse = '; ')))

  # Available analyses
  with_peaks = all(!is.null(x$peaks_analysis))
  with_CCF = all(!is.null(x$CCF_estimates))
  with_smoothing = all(!is.null(x$before_smoothing))
  with_arm_frag = all(!is.null(x$arm_fragmentation))
  with_wg_frag = all(!is.null(x$wg_fragmentation))

  if(with_peaks | with_CCF | with_smoothing | with_arm_frag | with_wg_frag) cat('\n')

  if(with_peaks)
  {
    cli::cli_alert_success("QC via peak detection available, score: {.value {x$peaks_analysis$score}}.")
    pio::pioDisp(x$peaks_analysis$matches)
  }

  if(with_CCF)
    cli::cli_alert_success("Cancer Cell Fraction (CCF) data available for karyotypes: {.value {names(x$CCF_estimates)}}.")

  if(with_smoothing)
    cli::cli_alert_success("These segments are smoothed; before smoothing there were {.value {x$before_smoothing$n_cna}} segments.")

  if(with_arm_frag)
    cli::cli_alert_success("Arm-level fragmentation analysis: {.value {sum(x$arm_fragmentation$table$significant)}} segments overfragmented.")

  if(with_wg_frag)
    {
    cond = x$wg_fragmentation$is_overfragmented
    p = round(x$wg_fragmentation$pvalue, 5)
    cli::cli_alert_success("Whole-genome fragmentation analysis: p = {.value {p}}: {.value {ifelse(cond, crayon::red('overfragmented'), crayon::green('not overfragmented'))}}.")
    }
}
