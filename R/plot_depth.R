#' Plot a genome-wide scatter plot of coverage.
#'
#'  @description Plot a genome-wide scatter plot of mutation depths, downsampled
#'  if required (annotates the used proportion).
#'
#' @param x An object of class \code{cnaqc}, created by the \code{init} function.
#' @param N Mutations to use, randomly sampled.
#' @param chromosomes The chromosome to use for this plot.
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna,example_dataset_CNAqc$purity)
#'
#' plot_depth(x)
#'
#' plot_depth(x, N = 100)
#' plot_depth(x, N = 1000)
#' plot_depth(x, N = 10000)
plot_depth = function(x, N = 5000, chromosomes = paste0('chr', c(1:22, 'X', 'Y')))
{
  stopifnot(inherits(x, 'cnaqc'))

  mutations = x$snvs %>%
    filter(chr %in% chromosomes) %>%
    relative_to_absolute_coordinates

  # X-range
  data('chr_coordinates_hg19', package = 'CNAqc')
  low = min(chr_coordinates_hg19$from)
  upp = max(chr_coordinates_hg19$to)

  med_DP = median(mutations$DP)

  quant = quantile(mutations$DP, probs = c(.1, .99))

  N_all = nrow(mutations)
  mutations = mutations %>%
    filter(DP > quant[1], DP < quant[2]) %>%
    sample_n(N)

  maxY = max(mutations$DP) * .9
  label_maxY = paste0("N = ", N, ' (', round(N/N_all * 100), '%)')

  ggplot(mutations,
                  aes(x = from, y = DP)) +
    scale_fill_viridis_c() +
    geom_point(size = .05) +
    my_ggplot_theme() +
    xlim(low, upp) +
    # ylim()
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    # ggpubr::rotate_y_text() +
    labs(y = "DP") +
    geom_hline(yintercept = med_DP, size = .4, linetype = 'dashed', color = 'darkred') +
    guides(fill = FALSE) +
    annotate("label", fill = 'white', x = upp, y = maxY, label = label_maxY, size = 2, hjust = 1)
}
