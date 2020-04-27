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
#' plot_gw_counts(x)
plot_gw_counts = function(x, chromosomes = paste0('chr', c(1:22, 'X', 'Y')))
{
  stopifnot(inherits(x, 'cnaqc'))

  mutations = CNAqc:::relative_to_absolute_coordinates(
    x,
    x$snvs %>% dplyr::filter(chr %in% chromosomes))

  bl_plot = CNAqc:::blank_genome(ref = x$reference_genome, chromosomes = chromosomes) +
    ylim(NA, NA)

  # # X-range
  # reference_genome = CNAqc:::get_reference(x$reference_genome)
  # low = min(reference_genome$from)
  # upp = max(reference_genome$to)

  # Histogram of mutation counts with 1 megabase bins
  binsize = 1e6

  hplot = bl_plot +
    geom_histogram(data = mutations, aes(x = from, y = ..count..), binwidth = binsize, fill = 'black') +
    CNAqc:::my_ggplot_theme() +
    theme(
      # axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    labs(y = 'n')

  hplot

}

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
#' plot_gw_depth(x)
#'
#' plot_gw_depth(x, N = 100)
#' plot_gw_depth(x, N = 1000)
#' plot_gw_depth(x, N = 10000)
plot_gw_depth = function(x, N = 5000, chromosomes = paste0('chr', c(1:22, 'X', 'Y')))
{
  stopifnot(inherits(x, 'cnaqc'))

  mutations = CNAqc:::relative_to_absolute_coordinates(
    x,
    x$snvs %>% dplyr::filter(chr %in% chromosomes))

  bl_plot = CNAqc:::blank_genome(ref = x$reference_genome, chromosomes = chromosomes) +
    ylim(NA, NA)

  # X-range
  # reference_genome = CNAqc:::get_reference(x$reference_genome)
  # low = min(reference_genome$from)
  # upp = max(reference_genome$to)

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

  dp = bl_plot +
    geom_point(data = mutations,
               aes(x = from, y = DP),
               size = .05 * cex_opt) +
    scale_fill_viridis_c() +
    # xlim(low, upp) +
    CNAqc:::my_ggplot_theme() +
    theme(
      # axis.text.x = element_blank(),
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
#' plot_gw_vaf(x)
#'
#' plot_gw_vaf(x, N = 100)
#' plot_gw_vaf(x, N = 1000)
#' plot_gw_vaf(x, N = 10000)
plot_gw_vaf = function(x, N = 5000, chromosomes = paste0('chr', c(1:22, 'X', 'Y')))
{
  stopifnot(inherits(x, 'cnaqc'))

  mutations = CNAqc:::relative_to_absolute_coordinates(
    x,
    x$snvs %>% dplyr::filter(chr %in% chromosomes))

  bl_plot = CNAqc:::blank_genome(ref = x$reference_genome, chromosomes = chromosomes) +
    ylim(NA, NA)

  # # X-range
  # reference_genome = CNAqc:::get_reference(x$reference_genome)
  # low = min(reference_genome$from)
  # upp = max(reference_genome$to)

  # VAF stats
  med_VAF = median(mutations$VAF)
  quant = quantile(mutations$DP, probs = c(.1, .99))

  N_all = nrow(mutations)
  mutations = mutations %>%
    filter(DP > quant[1], DP < quant[2])

  if(nrow(mutations) > N)
    mutations = mutations %>% sample_n(N)

  maxY = max(mutations$VAF) * .9
  label_maxY = paste0("N = ", N, ' (', round(N/N_all * 100), '%)')

  vaf = bl_plot +
    geom_point(data = mutations, aes(x = from, y = VAF), size = .05) +
    CNAqc:::my_ggplot_theme() +
    # scale_x_continuous(
    #   breaks = c(0, upp),
    #   labels = c("", "")
    # ) +
    # xlim(low, upp) +
    theme(
      # axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    ylim(0,1) +
    labs(y = "VAF") +
    geom_hline(yintercept = med_VAF, size = .4, linetype = 'dashed', color = 'darkred') +
    guides(fill = FALSE)


  # Simulate an internal legend
  L = ggplot_build(vaf)$layout$panel_params[[1]]
  Lx = abs(L$x.range[2] - L$x.range[1]) * .85

  vaf + annotate("label", fill = 'white', x = Lx, y = maxY, label = label_maxY, size = 2, hjust = 1)
}

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
#' # Compute CCF
#' x = compute_CCF(x)
#' plot_gw_ccf(x)
#'
#' plot_gw_ccf(x, N = 100)
#' plot_gw_ccf(x, N = 1000)
#' plot_gw_ccf(x, N = 10000)
plot_gw_ccf = function(x, N = 5000, chromosomes = paste0('chr', c(1:22, 'X', 'Y')))
{
  stopifnot(inherits(x, 'cnaqc'))

  with_CCF = all(!is.null(x$CCF_estimates))

  if(!with_CCF) {
    warning("Input does not have CCF estimates, see ?compute_CCF to determine CCF values.")
    return(CNAqc:::eplot())
  }

  mutations = CNAqc:::relative_to_absolute_coordinates(
    x,
    CNAqc::CCF(x) %>% dplyr::filter(chr %in% chromosomes)
    )

  bl_plot = CNAqc:::blank_genome(ref = x$reference_genome, chromosomes = chromosomes) +
    ylim(NA, NA)

  # # X-range
  # reference_genome = CNAqc:::get_reference(x$reference_genome)
  # low = min(reference_genome$from)
  # upp = max(reference_genome$to)

  # CCF stats
  med_CCF = median(mutations$CCF)
  quant = quantile(mutations$DP, probs = c(.1, .99))

  N_all = nrow(mutations)
  mutations = mutations %>%
    filter(DP > quant[1], DP < quant[2])

  if(nrow(mutations) > N)
    mutations = mutations %>% sample_n(N)

  maxY = max(mutations$VAF) * .9
  label_maxY = paste0("N = ", N, ' (', round(N/N_all * 100), '%)')

  ccf = ggplot(mutations,
               aes(x = from, y = CCF)) +
    geom_point(size = .05) +
    CNAqc:::my_ggplot_theme() +
    # scale_x_continuous(
    #   breaks = c(0, upp),
    #   labels = c("", "")
    # ) +
    # xlim(low, upp) +
    theme(
      # axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    labs(y = "CCF") +
    geom_hline(yintercept = med_CCF, size = .4, linetype = 'dashed', color = 'darkred') +
    guides(fill = FALSE)


  # Simulate an internal legend
  L = ggplot_build(ccf)$layout$panel_params[[1]]
  Lx = abs(L$x.range[2] - L$x.range[1]) * .85

  ccf + annotate("label", fill = 'white', x = Lx, y = maxY, label = label_maxY, size = 2, hjust = 1)
}


