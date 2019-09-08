#' Title
#'
#' @param x
#' @param digits
#'
#' @return
#' @export
#'
#' @examples
plot_peaks_analysis = function(x, digits = 4)
{
  stopifnot(inherits(x, "cnaqc"))

  with_peaks = all(!is.null(x$peaks_analysis))
  if(!with_peaks){
    stop("Input does not have peaks, see ?peaks_analysis to run peaks analysis.")
  }

  detections = x$peaks_analysis$plots

    ggpubr::ggarrange(
      plotlist = detections,
      nrow = 1,
      ncol = length(detections)
    )


}
