#' Title
#'
#' @param x
#' @param N
#' @param chromosomes
#'
#' @return
#' @export
#'
#' @examples
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
    filter(DP > quant[1], DP < quant[2]) %>%
    sample_n(N)

  maxY = max(mutations$VAF) * .9
  label_maxY = paste0("N = ", N, ' (', round(N/N_all * 100), '%)')

  vaf = ggplot(mutations,
         aes(x = from, y = VAF)) +
    geom_point(size = .05) +
    my_ggplot_theme() +
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
