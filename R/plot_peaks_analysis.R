#' Plot the results of peak analysis.
#'
#'  @description Results from \code{analyze_peaks} can be visualised with this
#'  function, which arranges the plots of each karyotype in a figure via \code{ggpubr}.
#'  Each karyotype shows the data, the estimated density, the peaks (selected and
#'  discarded), and the fit with shaded matching area.
#'
#' @param x An object of class \code{cnaqc}, where function \code{analyze_peaks} has
#' been computed.
#'
#' @return A \code{ggpubr} object for an assembled figure.
#' @export
#'
#' @import ggpubr
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna,example_dataset_CNAqc$purity)
#'
#' x = analyze_peaks(x)
#' plot_peaks_analysis(x)
plot_peaks_analysis = function(x)
{
  stopifnot(inherits(x, "cnaqc"))

  with_peaks = all(!is.null(x$peaks_analysis))
  if(!with_peaks){
    stop("Input does not have peaks, see ?peaks_analysis to run peaks analysis.")
  }

  detections = x$peaks_analysis$plots

  pl = suppressWarnings(ggpubr::ggarrange(
    plotlist = detections,
    nrow = 1,
    ncol = length(detections)
    ))

  qc = ifelse(x$peaks_analysis$QC == 'PASS', 'forestgreen', 'indianred3')

  pl = pl +
    theme(title = element_text(color = qc),
          panel.border = element_rect(
            colour = qc,
            fill = NA
          ))

  pl
}
