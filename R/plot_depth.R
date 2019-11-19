#' Plot a genome-wide scatter plot of coverage.
#'
#'  @description Plot a genome-wide scatter plot of mutation depths, downsampled
#'  if required (annotates the used proportion).
#'
#' @param x An object of class \code{cnaqc}, created by the \code{init} function.
#' @param N Mutations to use, randomly sampled.
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
#' plot_depth(x)
#'
#' plot_depth(x, N = 100)
#' plot_depth(x, N = 1000)
#' plot_depth(x, N = 10000)
plot_depth = function(x, N = 5000, chromosomes = paste0('chr', c(1:22, 'X', 'Y')), annotate_chromosomes = FALSE)
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
    filter(DP > quant[1], DP < quant[2])

  if(nrow(mutations) > N)
    mutations = mutations %>% sample_n(N)

  maxY = round(max(mutations$DP)  / 20) * 20
  minY = floor(min(mutations$DP)  / 20) * 20
  label_maxY = paste0("N = ", N, ' (', round(N/N_all * 100), '%)')
  
  base_plot = ggplot()
  if (annotate_chromosomes == TRUE) {
    base_plot = blank_genome(chromosomes, y = minY)  
  }
  
  depthplot = base_plot +
    scale_fill_viridis_c() +
    geom_point(data = mutations, aes(x = from, y = DP), size = .05) +
    my_ggplot_theme() +
    xlim(low, upp) +
    ylim(c(minY, maxY)) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    labs(y = "DP") +
    geom_hline(yintercept = med_DP, size = .4, linetype = 'dashed', color = 'darkred') +
    guides(fill = FALSE) +
    annotate("label", fill = 'white', x = upp, y = maxY * 0.9, label = label_maxY, size = 2, hjust = 1)
  
  depthplot
}
