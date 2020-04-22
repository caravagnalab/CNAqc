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

  mutations = CNAqc:::relative_to_absolute_coordinates(
    x,
    x$snvs %>% dplyr::filter(chr %in% chromosomes))

  # X-range
  reference_genome = CNAqc:::get_reference(x$reference_genome)
  low = min(reference_genome$from)
  upp = max(reference_genome$to)

  med_DP = median(mutations$DP)

  quant = quantile(mutations$DP, probs = c(.1, .99))

  N_all = nrow(mutations)
  mutations = mutations %>%
    filter(DP > quant[1], DP < quant[2])

  if(nrow(mutations) > N)
    mutations = mutations %>% sample_n(N)

  maxY = max(mutations$DP) * .9
  label_maxY = paste0("N = ", N, ' (', round(N/N_all * 100), '%)')

  cex_opt = getOption('CNAqc_cex', default = 1)

  dp = ggplot(mutations,
                  aes(x = from, y = DP)) +
    scale_fill_viridis_c() +
    xlim(low, upp) +
    geom_point(size = .05 * cex_opt) +
    CNAqc:::my_ggplot_theme() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    labs(y = "DP") +
    geom_hline(yintercept = med_DP, size = .4, linetype = 'dashed', color = 'darkred') +
    guides(fill = FALSE)

  # Simulate an internal legend
  L = ggplot_build(dp)$layout$panel_params[[1]]
  Lx = abs(L$x.range[2] - L$x.range[1]) * .85

  dp + annotate("label", fill = 'white', x = Lx, y = maxY, label = label_maxY, size = 2, hjust = 1)
}
