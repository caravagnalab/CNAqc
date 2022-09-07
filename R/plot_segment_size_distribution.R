#' Plot the length of clonal simple CNAs.
#'
#' @description
#'
#' This is the histogram of the number of bases (length)
#' in each copy number segment, for clonal simple CNAs.
#'
#' @param x A CNAqc object.
#' @param annotate_cut A custom vertical line annotation, by default at `1e7` megabases.
#'
#' @return A \code{ggplot2} plot.
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(mutations = example_dataset_CNAqc$mutations, cna = example_dataset_CNAqc$cna, purity = example_dataset_CNAqc$purity)
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
