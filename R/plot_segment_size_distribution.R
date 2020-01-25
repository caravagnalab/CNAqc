#' Title
#'
#' @param x
#' @param annotate_cut
#'
#' @return
#' @export
#'
#' @examples
plot_segment_size_distribution = function(x, annotate_cut = 1e7)
{
  x$cna %>%
    ggplot() +
    geom_histogram(aes(length), bins = 30) +
    CNAqc:::my_ggplot_theme() +
    labs(
      x = "Segments size",
      title = "Segments size distribution",
      captionn = paste0("Dashed line: ", annotate_cut)
    ) +
    geom_vline(xintercept = annotate_cut, color = 'red', linetype = 'dashed', size = .3)
}
