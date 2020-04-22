#' Plot a genome-wide histogram of mutation counts.
#'
#'  @description Plot a genome-wide histogram of mutation counts, binned
#' every one megabase (10e6 positions according to the hg19 reference).
#'
#' @param x An object of class \code{cnaqc}, created by the \code{init} function.
#' @param chromosomes The chromosome to use for this plot.
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna,example_dataset_CNAqc$purity)
#'
#' plot_counts(x)
plot_counts = function(x, chromosomes = paste0('chr', c(1:22, 'X', 'Y')))
{
  stopifnot(inherits(x, 'cnaqc'))

  mutations = CNAqc:::relative_to_absolute_coordinates(
    x,
    x$snvs %>% dplyr::filter(chr %in% chromosomes))

  # X-range
  reference_genome = CNAqc:::get_reference(x$reference_genome)
  low = min(reference_genome$from)
  upp = max(reference_genome$to)

  # Histogram of mutation counts with 1 megabase bins
  binsize = 1e6

  hplot = ggplot(mutations, aes(x = from)) +
    geom_histogram(aes(y = ..count..), binwidth = binsize, fill = 'black') +
    my_ggplot_theme() +
    xlim(low, upp) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    labs(y = 'n')

  hplot

}
