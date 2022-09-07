#' Plot the segments' size distribution.
#'
#' @description
#'
#' This is the histogram of the number of bases (length)
#' in each copy number segment.
#'
#' @param x An object of class \code{cnaqc}, created by the \code{init} function.
#' @param annotate_cut A custom vertical line annotation.
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(example_dataset_CNAqc$mutations, example_dataset_CNAqc$cna,example_dataset_CNAqc$purity)
#'
#' plot_segment_size_distribution(x)
plot_segment_size_distribution = function(x, annotate_cut = 1e7)
{
  x$cna %>%
    ggplot2::ggplot() +
    ggplot2::geom_histogram(ggplot2::aes(length), bins = 30) +
    CNAqc:::my_ggplot_theme() +
    ggplot2::labs(
      x = "Segments size",
      title = "Segments size distribution",
      caption = paste0("Dashed line: ", annotate_cut)
    ) +
    ggplot2::geom_vline(xintercept = annotate_cut, color = 'red', linetype = 'dashed', size = .3)
}
