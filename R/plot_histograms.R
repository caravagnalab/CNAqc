#' Title
#'
#' @param x
#' @param which
#' @param karyotypes
#'
#' @return
#' @export
#'
#' @examples
plot_data_histogram = function(x, which = 'VAF', karyotypes = c("1:0", '1:1', '2:0', '2:1', '2:2'))
{
  stopifnot(inherits(x, 'cnaqc'))

  if (which == 'VAF')
    return(plot_VAF_data(x, karyotypes = karyotypes))
  if (which == 'DP')
    return(plot_DP_data(x, karyotypes = karyotypes))
  if (which == 'CCF')
    return(plot_CCF_data(x, karyotypes = karyotypes))

  stop("Which one of VAF, DP or CCF.")
}


# Plot a histogram of CCF data
plot_CCF_data = function(x, karyotypes = c("1:0", '1:1', '2:0', '2:1', '2:2'))
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

  ggplot(data = ccf_data,
         aes(CCF, fill = karyotype)) +
    geom_histogram(binwidth = 0.01) +
    xlim(0, max(ccf_data$CCF, na.rm = T) %>% ceiling) +
    CNAqc:::my_ggplot_theme() +
    labs(title = "CCF values",
         caption = paste0("n = ", nrow(ccf_data))) +
    scale_fill_manual(values = CNAqc:::get_karyotypes_colors(unique(ccf_data$karyotype)))
}


# Plot a histogram of VAF data
plot_VAF_data = function(x, karyotypes = c("1:0", '1:1', '2:0', '2:1', '2:2'))
{
  stopifnot(inherits(x, 'cnaqc'))

  karyotypes = c("1:0", '1:1', '2:0', '2:1', '2:2')

  # VAF
  raw_muts = x$snvs %>%
    dplyr::mutate(karyotype = ifelse(karyotype %in% karyotypes, karyotype, "other"))

  ggplot(data = raw_muts, aes(VAF, fill = karyotype)) +
    geom_histogram(binwidth = 0.01) +
    xlim(0, 1) +
    CNAqc:::my_ggplot_theme() +
    labs(title = "Raw VAF",
         caption = paste0(
           "n = ",
           nrow(raw_muts),
           "; VAF < 0.05 (5%) = ",
           sum(x$snvs$VAF < 0.05)
         )) +
    scale_fill_manual(values = CNAqc:::get_karyotypes_colors(unique(raw_muts$karyotype)))


}

# Plot the same for DP
plot_DP_data = function(x, karyotypes = c("1:0", '1:1', '2:0', '2:1', '2:2'))
{
  stopifnot(inherits(x, 'cnaqc'))

  karyotypes = c("1:0", '1:1', '2:0', '2:1', '2:2')

  # VAF
  raw_muts = x$snvs %>%
    dplyr::mutate(karyotype = ifelse(karyotype %in% karyotypes, karyotype, "other"))

  ggplot(data = raw_muts, aes(DP, fill = karyotype)) +
    geom_histogram(bins = 100) +
    scale_x_log10() +
    CNAqc:::my_ggplot_theme() +
    labs(title = "Sequencing depth (coverage)",
         caption = paste0("Median coverage ", median(raw_muts$DP), 'x')) +
    geom_vline(
      xintercept = median(raw_muts$DP),
      color = 'orange',
      linetype = 'dashed',
      size = .3
    ) +
    scale_fill_manual(values = CNAqc:::get_karyotypes_colors(unique(raw_muts$karyotype)))
}
