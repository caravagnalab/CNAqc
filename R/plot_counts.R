#' Plot a genome-wide histogram of mutation counts.
#'
#'  @description Plot a genome-wide histogram of mutation counts, binned
#' every one megabase (10e6 positions according to the hg19 reference).
#'
#' @param x An object of class \code{cnaqc}, created by the \code{init} function.
#' @param chromosomes The chromosome to use for this plot.
#' @param annotate_chromosomes Boolean value specifying if chromosome should be annotated or not, default = FALSE
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna,example_dataset_CNAqc$purity)
#'
#' plot_counts(x)
plot_counts = function(x, chromosomes = paste0('chr', c(1:22, 'X', 'Y')), annotate_chromosomes = FALSE)
{
  stopifnot(inherits(x, 'cnaqc'))

  mutations = x$snvs %>%
    filter(chr %in% chromosomes) %>%
    relative_to_absolute_coordinates

  # X-range
  data('chr_coordinates_hg19', package = 'CNAqc')
  low = min(chr_coordinates_hg19$from)
  upp = max(chr_coordinates_hg19$to)

  # Histogram of mutation counts with 1 megabase bins
  binsize = 1e6
  
  base_plot = ggplot()
  if (annotate_chromosomes == TRUE) {
    base_plot = blank_genome(chromosomes, y = -1, chrlinecolour = "firebrick4")  
  }

  hplot = base_plot +
    geom_histogram(data = mutations, aes(x = from, y = ..count..), binwidth = binsize, fill = 'black') +
    my_ggplot_theme() +
    xlim(low, upp) +
    # scale_x_continuous(low, upp - upp*binsize) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    # ggpubr::rotate_y_text() +
    # scale_y_continuous(breaks = scales::pretty_breaks(n = 2),
    #                    limits = c(0, m)) +
    labs(y = 'n')

  hplot

}
