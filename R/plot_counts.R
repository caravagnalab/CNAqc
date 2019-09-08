#' Title
#'
#' @param x
#' @param chromosomes
#'
#' @return
#' @export
#'
#' @examples
plot_counts = function(x, chromosomes = paste0('chr', c(1:22, 'X', 'Y')))
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

  hplot = ggplot(mutations, aes(x = from)) +
    geom_histogram(aes(y = ..count..), binwidth = binsize, fill = 'black') +
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
