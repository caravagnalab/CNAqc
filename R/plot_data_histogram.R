#' Plot the data histogram.
#'
#' @description
#'
#' This function plots the histogram of any of the following:
#' the Variant Allele Frequency (VAF), the depth of sequencing
#' (DP), the nuber of reads with the variant (NV) and the
#' Cancer Cell Fractions (CCF) esimates (that must be computed
#' before.
#'
#' The plot can be subset by selecting only data in certain segments,
#' identified by the karyotype shortname, and will be coloured accordingly.
#'
#' @param x An object of class \code{cnaqc}, created by the \code{init} function.
#' @param which One of \code{"VAF"}, \code{"DP"}, \code{"NV"} or \code{"CCF"}.
#' @param karyotype A list of karyotype ids in \code{"Major:minor"} notation
#' (e.g., \code{"1:1", "2,1", ...}) that will be retained ofr the plot.
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @examples
#' #' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna,example_dataset_CNAqc$purity)
#' plot_data_histogram(x)
#' plot_data_histogram(x, which = 'DP')
#' plot_data_histogram(x, which = 'DP', karyotypes = '2:2')
#'
#' x = compute_CCF(x)
#' plot_data_histogram(x, which = 'CCF')
plot_data_histogram = function(x,
                               which = 'VAF',
                               karyotypes = c("1:0", '1:1', '2:0', '2:1', '2:2'))
{
  stopifnot(inherits(x, 'cnaqc'))

  if(!(which %in% c("NV", "DP", "VAF", "CCF")))
    stop("'which' must be one of 'VAF', 'DP', 'NV' or 'CCF'.")

  # By case
  plot_f = CNAqc:::eplot()

  if (which == 'VAF') plot_f = CNAqc:::plot_VAF_data(x, karyotypes = karyotypes) + 
    facet_wrap(~type, ncol = 1, scales = 'free_y')
  
  if (which == 'DP') plot_f = CNAqc:::plot_DP_data(x, karyotypes = karyotypes) +
    facet_wrap(~type, ncol = 1, scales = 'free_y')
  
  if (which == 'NV') plot_f = CNAqc:::plot_NV_data(x, karyotypes = karyotypes) +
    facet_wrap(~type, ncol = 1, scales = 'free_y')
  
  with_CCF = all(!is.null(x$CCF_estimates))
  if (which == 'CCF') {
    plot_f = CNAqc:::plot_CCF_data(x, karyotypes = karyotypes) 
    
    if(with_CCF) plot_f = plot_f + facet_wrap(~type, ncol = 1, scales = 'free_y')
  }

  if(!CNAqc:::has_driver_data(x)) return(plot_f)

  # Drivers
  if(
    which != "CCF" | with_CCF
  )
  plot_f = CNAqc:::annotate_drivers_to_histogram(
    drivers_list = get_drivers(x,  which = ifelse(which %in% c("VAF", "CCF"), which, 'VAF')) %>%
      dplyr::mutate(karyotype = ifelse(
        karyotype %in% karyotypes, karyotype, "other"
      )),
    p = plot_f,
    which = which
  )

  return(plot_f)
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

  ggplot(data = ccf_data,
         aes(CCF, fill = karyotype)) +
    geom_histogram(binwidth = 0.01) +
    xlim(-0.01, 0.01 + max(ccf_data$CCF, na.rm = T) %>% ceiling) +
    CNAqc:::my_ggplot_theme() +
    labs(title = bquote("CCF (" * bold(.(meth)) * ')'),
         caption = paste0("n = ", nrow(ccf_data))) +
    scale_fill_manual(values = CNAqc:::get_karyotypes_colors(unique(ccf_data$karyotype)))
}


# Plot a histogram of VAF data
plot_VAF_data = function(x,
                         karyotypes = c("1:0", '1:1', '2:0', '2:1', '2:2'))
{
  stopifnot(inherits(x, 'cnaqc'))

  # VAF
  raw_muts = x$snvs %>%
    dplyr::mutate(karyotype = ifelse(karyotype %in% karyotypes, karyotype, "other"))

  ggplot(data = raw_muts, aes(VAF, fill = karyotype)) +
    geom_histogram(binwidth = 0.01) +
    xlim(-0.01, 1.01) +
    CNAqc:::my_ggplot_theme() +
    labs(title = "VAF",
         caption = paste0(
           "n = ",
           nrow(raw_muts),
           "; VAF < 0.05 (5%) = ",
           sum(x$snvs$VAF < 0.05)
         )) +
    scale_fill_manual(values = CNAqc:::get_karyotypes_colors(unique(raw_muts$karyotype))) 
}

# Plot the same for DP
plot_DP_data = function(x,
                        karyotypes = c("1:0", '1:1', '2:0', '2:1', '2:2'))
{
  stopifnot(inherits(x, 'cnaqc'))

  # VAF
  raw_muts = x$snvs %>%
    dplyr::mutate(karyotype = ifelse(karyotype %in% karyotypes, karyotype, "other"))

  ggplot(data = raw_muts, aes(DP, fill = karyotype)) +
    geom_histogram(bins = 100) +
    scale_x_log10() +
    CNAqc:::my_ggplot_theme() +
    labs(title = "Sequencing depth",
         caption = paste0("Median DP ", median(raw_muts$DP), 'x')) +
    geom_vline(
      xintercept = median(raw_muts$DP),
      color = 'black',
      linetype = 'dashed',
      size = .5
    ) +
    scale_fill_manual(values = CNAqc:::get_karyotypes_colors(unique(raw_muts$karyotype)))
}

# Plot the same for NV
plot_NV_data = function(x,
                        karyotypes = c("1:0", '1:1', '2:0', '2:1', '2:2'))
{
  stopifnot(inherits(x, 'cnaqc'))

  # VAF
  raw_muts = x$snvs %>%
    dplyr::mutate(karyotype = ifelse(karyotype %in% karyotypes, karyotype, "other"))

  ggplot(data = raw_muts, aes(NV, fill = karyotype)) +
    geom_histogram(bins = 100) +
    scale_x_log10() +
    CNAqc:::my_ggplot_theme() +
    labs(title = "Variant reads",
         caption = paste0("Median NV ", median(raw_muts$NV), 'x')) +
    geom_vline(
      xintercept = median(raw_muts$NV),
      color = 'black',
      linetype = 'dashed',
      size = .5
    ) +
    scale_fill_manual(values = CNAqc:::get_karyotypes_colors(unique(raw_muts$karyotype)))
}

