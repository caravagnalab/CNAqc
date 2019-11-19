#' Plot a genome-wide scatter plot of VAF
#'
#'  @description Plot a genome-wide scatter plot of mutation VAFs, downsampled
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
#' plot_vaf(x)
#'
#' plot_vaf(x, N = 100)
#' plot_vaf(x, N = 1000)
#' plot_vaf(x, N = 10000)
plot_vaf = function(x, N = 5000, chromosomes = paste0('chr', c(1:22, 'X', 'Y')))
{
  stopifnot(inherits(x, 'cnaqc'))

  mutations = x$snvs %>%
    filter(chr %in% chromosomes) %>%
    relative_to_absolute_coordinates

  # X-range
  data('chr_coordinates_hg19', package = 'CNAqc')
  low = min(chr_coordinates_hg19$from)
  upp = max(chr_coordinates_hg19$to)

  med_VAF = median(mutations$VAF)

  quant = quantile(mutations$DP, probs = c(.1, .99))

  N_all = nrow(mutations)
  mutations = mutations %>%
    filter(DP > quant[1], DP < quant[2])

  if(nrow(mutations) > N)
    mutations = mutations %>% sample_n(N)


  maxY = max(mutations$VAF) * .9
  label_maxY = paste0("N = ", N, ' (', round(N/N_all * 100), '%)')
  
  base_plot = blank_genome(chromosomes, y = 0.025)
  
  vaf = vaf = base_plot + 
    geom_point(data = mutations, aes(x = from, y = VAF), size = .05) +
    my_ggplot_theme()  +
    xlim(low, upp) +
    # ylim()
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    ylim(0,1) +
    # ggpubr::rotate_y_text() +
    labs(y = "VAF") +
    geom_hline(yintercept = med_VAF, size = .4, linetype = 'dashed', color = 'darkred') +
    guides(fill = FALSE) +
    annotate("label", fill = 'white', x = upp, y = maxY, label = label_maxY, size = 2, hjust = 1)

  vaf
}
