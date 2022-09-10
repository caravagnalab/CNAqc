#' Plot the read-counts data histograms.
#'
#' @description
#'
#' This function plots the histogram of any of the following:
#'
#' * the Variant Allele Frequency (VAF),
#' * the depth of sequencing (DP),
#' * the number of reads with the variant (NV)
#' * the Cancer Cell Fractions (CCF) estimates (if available)
#'
#' The plot can be subset by selecting only mutations mapping to certain
#' karyotypes.
#'
#' @param x A CNAqc object.
#' @param which One of \code{"VAF"}, \code{"DP"}, \code{"NV"} or \code{"CCF"}.
#' @param karyotype A list of karyotypes in \code{"Major:minor"} notation
#' (e.g., \code{"1:1", "2,1", ...}) for the plot. By default \code{c("1:0", '1:1', '2:0', '2:1', '2:2')}
#' are used.
#'
#' @return A \code{ggplot2} plot.
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(mutations = example_dataset_CNAqc$mutations, cna = example_dataset_CNAqc$cna, purity = example_dataset_CNAqc$purity)
#'
#' # Default plot
#' plot_data_histogram(x)
#'
#' # Customised plots
#' plot_data_histogram(x, which = 'DP')
#' plot_data_histogram(x, which = 'DP', karyotypes = '2:2')
#'
#' # CCF computation and plotting
#' x = compute_CCF(x)
#' plot_data_histogram(x, which = 'CCF')
plot_data_histogram = function(x,
                               which = 'VAF',
                               karyotypes = c("1:0", '1:1', '2:0', '2:1', '2:2'))
{
  stopifnot(inherits(x, 'cnaqc'))

  if (!(which %in% c("NV", "DP", "VAF", "CCF")))
    stop("'which' must be one of 'VAF', 'DP', 'NV' or 'CCF'.")

  # By case
  plot_f = CNAqc:::eplot()

  if (which == 'VAF')
    plot_f = plot_VAF_data(x, karyotypes = karyotypes)

  if (which == 'DP')
    plot_f = plot_DP_data(x, karyotypes = karyotypes)

  if (which == 'NV')
    plot_f = plot_NV_data(x, karyotypes = karyotypes)

  with_CCF = all(!is.null(x$CCF_estimates))
  if (which == 'CCF') {
    plot_f = plot_CCF_data(x, karyotypes = karyotypes)

    if (with_CCF)
      plot_f = plot_f + ggplot2::facet_wrap( ~ type, ncol = 1, scales = 'free_y')
  }

  plot_f = plot_f + ggplot2::guides(fill = ggplot2::guide_legend(''))

  if (!CNAqc:::has_driver_data(x))
    return(plot_f)

  # Drivers
  if (which != "CCF" | with_CCF)
    plot_f = annotate_drivers_to_histogram(
      drivers_list = get_drivers(x,  which = ifelse(which %in% c("VAF", "CCF"), which, 'VAF')) %>%
        dplyr::mutate(karyotype = ifelse(
          karyotype %in% karyotypes, karyotype, "other"
        )),
      p = plot_f,
      which = which
    )

  return(plot_f +)
}


# Plot a histogram of CCF data
plot_CCF_data = function(x,
                         karyotypes = c("1:0", '1:1', '2:0', '2:1', '2:2'))
{
  stopifnot(inherits(x, 'cnaqc'))

  if (all(is.null(x$CCF_estimates)))
  {
    warning("Input does not have CCF estimates, see ?compute_CCF to determine CCF values.")
    return(CNAqc:::eplot())
  }

  # CCFs
  ccf_data = Reduce(dplyr::bind_rows,
                    lapply(x$CCF_estimates, function(x)
                      x$mutations)) %>%
    dplyr::mutate(karyotype = ifelse(karyotype %in% karyotypes, karyotype, "other")) %>%
    dplyr::filter(!is.na(CCF))

  # Whatever is fit
  meth = x$CCF_estimates[[1]]$QC_table$method

  ggplot2::ggplot(data = ccf_data,
                  ggplot2::aes(CCF, fill = karyotype)) +
    ggplot2::geom_histogram(binwidth = 0.01) +
    ggplot2::xlim(-0.01, 0.01 + max(ccf_data$CCF, na.rm = T) %>% ceiling) +
    CNAqc:::my_ggplot_theme() +
    ggplot2::labs(title = bquote("CCF (" * bold(.(meth)) * ')'),
                  caption = paste0("n = ", nrow(ccf_data))) +
    ggplot2::scale_fill_manual(values = CNAqc:::get_karyotypes_colors(unique(ccf_data$karyotype)))
}


# Plot a histogram of VAF data
plot_VAF_data = function(x,
                         karyotypes = c("1:0", '1:1', '2:0', '2:1', '2:2'))
{
  stopifnot(inherits(x, 'cnaqc'))

  # VAF
  # raw_muts = x$mutations %>%
  #   dplyr::mutate(karyotype = ifelse(karyotype %in% karyotypes, karyotype, "other"))

  raw_muts = Mutations(x) %>%
    dplyr::mutate(karyotype = ifelse(karyotype %in% karyotypes, karyotype, "other"))

  colors = get_karyotypes_colors(unique(raw_muts$karyotype))
  colors['subclonal'] = ggplot2::alpha('purple4', .7)

  ggplot2::ggplot(data = raw_muts, ggplot2::aes(VAF, fill = karyotype)) +
    ggplot2::geom_histogram(binwidth = 0.01) +
    ggplot2::xlim(-0.01, 1.01) +
    CNAqc:::my_ggplot_theme() +
    ggplot2::labs(title = "VAF",
                  caption = paste0(
                    "n = ",
                    nrow(raw_muts),
                    "; VAF < 0.05 (5%) = ",
                    sum(x$mutations$VAF < 0.05)
                  )) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::facet_wrap(type ~ paste("CNA", cna), scales = 'free_y')
}

# Plot the same for DP
plot_DP_data = function(x,
                        karyotypes = c("1:0", '1:1', '2:0', '2:1', '2:2'))
{
  stopifnot(inherits(x, 'cnaqc'))

  raw_muts = Mutations(x) %>%
    dplyr::mutate(karyotype = ifelse(karyotype %in% karyotypes, karyotype, "other"))

  colors = get_karyotypes_colors(unique(raw_muts$karyotype))
  colors['subclonal'] = ggplot2::alpha('purple4', .7)

  ggplot2::ggplot(data = raw_muts, ggplot2::aes(DP, fill = karyotype)) +
    ggplot2::geom_histogram(bins = 100) +
    ggplot2::scale_x_log10() +
    CNAqc:::my_ggplot_theme() +
    ggplot2::labs(title = "Sequencing depth",
                  caption = paste0("Median DP ", median(raw_muts$DP), 'x')) +
    ggplot2::geom_vline(
      xintercept = median(raw_muts$DP),
      color = 'black',
      linetype = 'dashed',
      size = .5
    ) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::facet_wrap(type ~ paste("CNA", cna), scales = 'free_y')

}

# Plot the same for NV
plot_NV_data = function(x,
                        karyotypes = c("1:0", '1:1', '2:0', '2:1', '2:2'))
{
  stopifnot(inherits(x, 'cnaqc'))

  # VAF
  raw_muts = Mutations(x) %>%
    dplyr::mutate(
      karyotype = case_when(
        karyotype %in% karyotypes ~ karyotype,
        is.na(karyotype) & cna == "subclonal" ~ 'subclonal',
        TRUE ~ "other"
      )
    )

  colors = get_karyotypes_colors(unique(raw_muts$karyotype))
  colors['subclonal'] = ggplot2::alpha('purple4', .7)

  ggplot2::ggplot(data = raw_muts, ggplot2::aes(NV, fill = karyotype)) +
    ggplot2::geom_histogram(bins = 100) +
    ggplot2::scale_x_log10() +
    CNAqc:::my_ggplot_theme() +
    ggplot2::labs(title = "Variant reads",
                  caption = paste0("Median NV ", median(raw_muts$NV), 'x')) +
    ggplot2::geom_vline(
      xintercept = median(raw_muts$NV),
      color = 'black',
      linetype = 'dashed',
      size = .5
    ) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::facet_wrap(type ~ paste("CNA", cna), scales = 'free_y')
}
