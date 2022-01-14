#' # Private function for sublconal CNAs peak detection and analysis
#' subclonal_peak_detector_closest_hit_match = function(chr,
#'                                                      from,
#'                                                      to,
#'                                                      mutations,
#'                                                      CCF_1,
#'                                                      CCF_2,
#'                                                      karyotype_1,
#'                                                      karyotype_2,
#'                                                      purity,
#'                                                      epsilon = 0.02)
#' {
#'   aux_peak_calling_subclonal = function(mutations) {
#'     xy_peaks = den = NULL
#'
#'     # Smoothed Gaussian kernel for VAF
#'     y = mutations %>% dplyr::pull(VAF)
#'
#'     den = density(y,
#'                   kernel = 'gaussian',
#'                   adjust = kernel_adjust,
#'                   na.rm = T)
#'     in_range = den$x >= min(y, na.rm = T) &
#'       den$x <= max(y, na.rm = T)
#'
#'     # den = density(y, kernel = 'gaussian', adjust = 0.5)
#'     # plot(den)
#'
#'     input_peakdetection = matrix(cbind(x = den$x[in_range], y = den$y[in_range]), ncol = 2)
#'     colnames(input_peakdetection) = c('x', 'y')
#'
#'     # Test 5 parametrisations of peakPick neighlim
#'     pks = Reduce(dplyr::bind_rows,
#'                  lapply(1:5,
#'                         function(n) {
#'                           pk = peakPick::peakpick(mat = input_peakdetection, neighlim = n)
#'                           input_peakdetection[pk[, 2], , drop = FALSE] %>% as.data.frame()
#'                         })) %>%
#'       as_tibble() %>%
#'       dplyr::arrange(x) %>%
#'       dplyr::mutate(x = round(x, 2), y = round(y, 2)) %>%
#'       dplyr::distinct(x, .keep_all = TRUE)
#'
#'     # print(pks)
#'     hst = hist(snvs$VAF, breaks = seq(0, 1, 0.01), plot = F)$counts
#'     pks$counts_per_bin = hst[round(pks$x * 100)]
#'
#'     # Heuristic to remove low-density peaks
#'     pks = pks %>%
#'       dplyr::mutate(discarded = y <= max(pks$y) * (1 / 20),
#'                     from = 'KDE')
#'
#'     return(list(peaks = pks, density = den))
#'   }
#'
#'   # Segments
#'   # segment_1 = subclonal_calls$karyotype_1[i]
#'   # segment_2 = subclonal_calls$karyotype_2[i]
#'   #
#'   # # Algo subclonal - subclones proportions
#'   # CCF_1 = subclonal_calls$CCF_1[i] # \psi
#'   # CCF_2 = subclonal_calls$CCF_2[i] # 1 - \psi
#'
#'   # Subclones ploidies
#'   ploidy_1 = strsplit(karyotype_1, ':')[[1]] %>% as.numeric() %>% sum # p_1
#'   ploidy_2 = strsplit(karyotype_2, ':')[[1]] %>% as.numeric() %>% sum # p_2
#'
#'   # Multiplicities expected - m_i_j ~ multiplicity j in clone i
#'   m_1_1 = m_1_2 = 1
#'   m_2_1 = m_2_2 =  2
#'
#'   # Common denominator
#'   denominator =
#'     (1 - purity) * 2 +
#'     ploidy_1 * CCF_1 * purity +
#'     ploidy_2 * CCF_2 * purity
#'
#'   # Numerator subclones -  v_i_j ~ reads from peak at multiplicity j in clone i
#'   v_1_1 = m_1_1 * CCF_1 * purity
#'   v_2_1 = m_2_1 * CCF_1 * purity
#'
#'   v_1_2 = m_1_2 * CCF_2 * purity
#'   v_2_2 = m_2_2 * CCF_2 * purity
#'
#'   # VAF peaks subclones ~ peak at multiplicity j in clone i
#'   peak_1_1 = v_1_1 / denominator
#'   peak_2_1 = v_2_1 / denominator
#'
#'   peak_1_2 = v_1_2 / denominator
#'   peak_2_2 = v_2_2 / denominator
#'
#'   # Filter unused peaks
#'   if (karyotype_1 %in% c("1:0", "1:1"))
#'     peak_2_1 = NA
#'   if (karyotype_1 %in% c("1:0", "1:1"))
#'     peak_2_2 = NA
#'
#'   # Data plot
#'   segment_plot =
#'     mutations %>%
#'     ggplot()  +
#'     geom_histogram(aes(x = VAF), binwidth = 0.01, fill = 'gray') +
#'     xlim(0, 1) +
#'     CNAqc:::my_ggplot_theme() +
#'     labs(
#'       title = paste(chr, from, to, sep = ':'),
#'       subtitle = paste0(
#'         karyotype_1,
#'         ' (',
#'         CCF_1 %>% round(2),
#'         ') - ',
#'         karyotype_2,
#'         ' (',
#'         CCF_2 %>% round(2),
#'         ') '
#'       )
#'     )
#'
#'   # Data peaks
#'   data_peaks = mutations %>%
#'     aux_peak_calling_subclonal()
#'
#'   # Check matching clone 1
#'   expected_1 =
#'     data.frame(peak = c(peak_1_1, peak_2_1),
#'                multiplicity = c(1, 2))
#'
#'   for (i in 1:nrow(expected_1))
#'   {
#'     expected_1$matched[i] = any(abs(data_peaks$peaks$x - expected_1$peak[i]) < epsilon)
#'   }
#'
#'   # Check matching clone 1
#'   expected_2 =
#'     data.frame(peak = c(peak_1_2, peak_2_2),
#'                multiplicity = c(1, 2))
#'
#'   for (i in 1:nrow(expected_2))
#'   {
#'     expected_2$matched[i] = any(abs(data_peaks$peaks$x - expected_2$peak[i]) < epsilon)
#'   }
#'
#'   # Overlay VAF density
#'   segment_plot = segment_plot +
#'     geom_line(
#'       data = data.frame(x = data_peaks$density$x,
#'                         y = data_peaks$density$y),
#'       aes(x = x, y = y),
#'       inherit.aes = FALSE
#'     )
#'
#'   # Plot peaks as points with overlaid shadow for matching
#'   for (i in 1:nrow(data_peaks$peaks))
#'     segment_plot = segment_plot +
#'     geom_rect(
#'       data = data.frame(
#'         xmin = data_peaks$peaks$x[i] - epsilon,
#'         xmax = data_peaks$peaks$x[i] + epsilon,
#'         ymin = 0,
#'         ymax = Inf
#'       ),
#'       aes(
#'         xmin = xmin,
#'         xmax = xmax,
#'         ymin = ymin,
#'         ymax = ymax
#'       ),
#'       inherit.aes = FALSE,
#'       fill = 'darkgoldenrod1',
#'       alpha = .3
#'     )
#'
#'   segment_plot = segment_plot +
#'     geom_point(
#'       data = data_peaks$peaks,
#'       aes(x = x, y = y),
#'       inherit.aes = FALSE,
#'       color = 'black'
#'     )
#'
#'   # Expected peaks clone 1
#'   if (any(!is.na(expected_1$peak)))
#'     segment_plot = segment_plot +
#'     geom_vline(
#'       xintercept = c(peak_1_1, peak_2_1),
#'       linetype = ifelse(expected_1$matched, 'solid', 'dashed'),
#'       color = 'indianred3',
#'       size = ifelse(expected_1$matched, 1, .5)
#'     )
#'
#'   # Expected peaks clone 2
#'   if (any(!is.na(expected_2$peak)))
#'     segment_plot = segment_plot +
#'     geom_vline(
#'       xintercept = c(peak_1_2, peak_2_2),
#'       linetype = ifelse(expected_2$matched, 'solid', 'dashed'),
#'       color = 'steelblue',
#'       size = ifelse(expected_2$matched, 1, .5)
#'     )
#'
#'   # Return tables
#'   expected_1$clone = 'Clone 1'
#'   expected_1$CCF = CCF_1
#'   expected_1$karyotype = karyotype_1
#'
#'   expected_2$clone = 'Clone 2'
#'   expected_2$CCF = CCF_2
#'   expected_2$karyotype = karyotype_2
#'
#'   list(
#'     table = dplyr::bind_rows(expected_1, expected_2) %>%
#'       dplyr::mutate(chr = chr, from = from, to = to) %>%
#'       dplyr::select(chr, from, to, everything()),
#'     plot = segment_plot
#'   ) %>%
#'     return()
#' }
#'
#'
#' # j = 7
#' # subclonal_peak_detector_closest_hit_match(
#' #   mutations = subclonal_calls$x[[j]],
#' #   CCF_1 = subclonal_calls$CCF_1[j],
#' #   CCF_2 = subclonal_calls$CCF_2[j],
#' #   karyotype_1 = subclonal_calls$karyotype_1[j],
#' #   karyotype_2 = subclonal_calls$karyotype_2[j],
#' #   chr = subclonal_calls$chr[j],
#' #   from = subclonal_calls$from[j],
#' #   to = subclonal_calls$to[j],
#' #   purity = purity,
#' #   epsilon = 0.03
#' # )
#'
#' # x = init(CNAqc::example_dataset_CNAqc$snvs, CNAqc::example_dataset_CNAqc$cna[-c(1:3), ], CNAqc::example_dataset_CNAqc$purity)
#' # subclonal_CNAs = CNAqc::example_dataset_CNAqc$cna[c(1:3), ]
#' # subclonal_CNAs$CCF_1 = runif(3)
#' # subclonal_CNAs$CCF_2 = 1 - subclonal_CNAs$CCF_1
#' # subclonal_CNAs$Major_1 = 2
#' # subclonal_CNAs$minor_1 = 1
#' # subclonal_CNAs$Major_2 = 2
#' # subclonal_CNAs$minor_2 = 2
#'
#' analyse_subclonal_peaks = function(x, mutations, subclonal_CNAs, epsilon)
#' {
#'   stopifnot(inherits(x, 'cnaqc'))
#'
#'   # Format subclonal_CNAs
#'   req_names = c("chr",
#'                 "from",
#'                 "to",
#'                 "CCF_1",
#'                 "CCF_2",
#'                 "Major_1",
#'                 "minor_1",
#'                 "Major_2",
#'                 "minor_2")
#'   missing_names = setdiff(req_names, colnames(subclonal_CNAs))
#'
#'   if (length(missing_names) > 0)
#'   {
#'     cli::cli_alert_danger(
#'       "The following columns are missing from your subclonal CNAs table: {.field {missing_names}}."
#'     )
#'     return(x)
#'   }
#'
#'   # Format mutations
#'   req_names = c("chr", "from", "to", "VAF")
#'   missing_names = setdiff(req_names, colnames(mutations))
#'
#'   if (length(missing_names) > 0)
#'   {
#'     cli::cli_alert_danger(
#'       "The following columns are missing from your mutations table: {.field {missing_names}}."
#'     )
#'     return(x)
#'   }
#'
#'   # Map all input - retain only supported peaks
#'   cli::cli_h1("Peak analysis for subclonal CNAs")
#'
#'   supported_peaks = c("1:0", "1:1", "2:0", "2:1", "2:2")
#'
#'   subclonal_CNAs = subclonal_CNAs %>%
#'     mutate(
#'       karyotype_1 = paste(Major_1, minor_1, sep = ':'),
#'       karyotype_2 = paste(Major_2, minor_2, sep = ':'),
#'       analysed =
#'         (karyotype_1 %in% supported_peaks) &
#'         (karyotype_2 %in% supported_peaks)
#'     )
#'
#'   if (any(!subclonal_CNAs$analysed))
#'   {
#'     cli::cli_alert_warning(
#'       "The following segments have subclonal CNAs that are not supported by CNAqc and will be rejected."
#'     )
#'     subclonal_CNAs %>% dplyr::filter(!analysed) %>% print()
#'     subclonal_CNAs = subclonal_CNAs %>% dplyr::filter(analysed)
#'   }
#'
#'   # Map mutations to the actual segments
#'   subclonal_CNAs$mutations = NULL
#'   for (i in 1:nrow(subclonal_CNAs))
#'   {
#'     subclonal_CNAs$mutations[i] = mutations %>%
#'       dplyr::filter(chr == subclonal_CNAs$chr[i],
#'                     from  >= subclonal_CNAs$from[i],
#'                     to <= subclonal_CNAs$to[i]) %>% list()
#'   }
#'
#'   if (any(sapply(subclonal_CNAs$mutations, nrow) >= 1))
#'   {
#'     to_remove = subclonal_CNAs[sapply(subclonal_CNAs$mutations, nrow) >= 1,]
#'     cli::cli_alert_warning("The following subclonal CNAs have <=1 mutations mapped and will not be analysed")
#'     to_remove %>% print()
#'
#'     subclonal_CNAs = subclonal_CNAs[sapply(subclonal_CNAs$mutations, nrow) >
#'                                       1,]
#'   }
#'
#'   # Actual analysis
#'   all_fits_plots = all_fits_table =  NULL
#'
#'   for (j in 1:nrow(subclonal_CNAs))
#'   {
#'     this_fit = subclonal_peak_detector_closest_hit_match(
#'       mutations = subclonal_CNAs$mutations[[j]],
#'       CCF_1 = subclonal_CNAs$CCF_1[j],
#'       CCF_2 = subclonal_CNAs$CCF_2[j],
#'       karyotype_1 = subclonal_CNAs$karyotype_1[j],
#'       karyotype_2 = subclonal_CNAs$karyotype_2[j],
#'       chr = subclonal_CNAs$chr[j],
#'       from = subclonal_CNAs$from[j],
#'       to = subclonal_CNAs$to[j],
#'       purity = x$purity,
#'       epsilon = 0.03
#'     )
#'
#'     all_fits_plots = append(all_fits_plots, list(this_fit$plot))
#'     all_fits_table = append(all_fits_table, list(this_fit$table))
#'   }
#'
#'   # Assembly plots
#'   all_fits_table = Reduce(bind_rows, all_fits_table) %>%
#'     mutate(id = paste(chr, from, to, sep = ":"))
#'
#'   tab_heatmap = ggplot(all_fits_table) +
#'     geom_tile(aes(y = id,
#'                   x = multiplicity,
#'                   fill = matched),
#'               width = .8,
#'               height = .8) +
#'     scale_fill_manual(values = c(
#'       `TRUE` = 'forestgreen',
#'       `FALSE` = 'indianred',
#'       `NA` = 'gray'
#'     )) +
#'     facet_wrap( ~ clone) +
#'     CNAqc:::my_ggplot_theme() +
#'     scale_y_discrete(limits = all_fits_table$id %>% unique %>% gtools::mixedsort(decreasing = TRUE)) +
#'     scale_x_discrete(limits = c(1:2) %>% paste) +
#'     guides(fill = guide_legend(ncol = 1))
#'
#'   karyo_colors = CNAqc:::get_karyotypes_colors(all_fits_table$karyotype %>% unique)
#'   karyo_colors = karyo_colors[all_fits_table$karyotype %>% unique]
#'
#'   karyo_heatmap = ggplot(all_fits_table) +
#'     geom_tile(aes(y = id,
#'                   x = clone,
#'                   fill = karyotype),
#'               width = .8,
#'               height = .8) +
#'     scale_fill_manual(values = karyo_colors) +
#'     CNAqc:::my_ggplot_theme() +
#'     scale_y_discrete(limits = all_fits_table$id %>% unique %>% gtools::mixedsort(decreasing = TRUE)) +
#'     guides(fill = guide_legend(ncol = 2)) +
#'     theme(axis.text.y = element_blank()) +
#'     labs(y = NULL)
#'
#'
#'   size_barplot = all_fits_table %>%
#'     ggplot() +
#'     geom_bar(aes(x = id, fill = clone, CCF),
#'              stat = 'identity',
#'              width = .8) +
#'     coord_flip()  +
#'     scale_x_discrete(limits = all_fits_table$id %>% unique %>% gtools::mixedsort(decreasing = TRUE)) +
#'     CNAqc:::my_ggplot_theme() +
#'     theme(axis.text.y = element_blank()) +
#'     labs(x = NULL, y = "Proportion") +
#'     guides(fill = guide_legend(ncol = 1))
#'
#'   figure = cowplot::plot_grid(
#'     tab_heatmap,
#'     karyo_heatmap,
#'     size_barplot,
#'     rel_widths = c(1, .2, .3),
#'     axis = 'tb',
#'     align = 'h',
#'     ncol = 3,
#'     nrow = 1
#'   )
#'
#'   x$subclonal = list(figure = list(figure = figure, fits = all_fits_plots),
#'                      calls = subclonal_CNAs,)
#'
#'   return(x)
#' }
#'
#' #' Parse Battenberg clonal/ subclonal calls
#' #'
#' #' @description Returns two tibbles in CNAqc-ready format, one with clonal calls
#' #' and the other with subclonal calls extracted by parsing the Battenberg format.
#' #'
#' #' @param x The input dataframe, required to have the following columns: \code{chr},
#' #' \code{from}, \code{to}, \code{battenberg_nMaj1_A}, \code{battenberg_nMin1_A},
#' #' \code{battenberg_nMaj2_A}, \code{battenberg_nMin2_A}, \code{battenberg_frac1_A},
#' #' and \code{battenberg_frac2_A}.
#' #'
#' #' @return Two tibbles in CNAqc-ready format.
#' #' @export
#' #'
#' #' @examples
#' #' \dontrun{
#' #' # Load some CSV results from Battenberg
#' #' x = read.csv(....) %>% parse_Battenberg()
#' #'
#' #' # Work with clonal calls (omitting mutations and other parameters here)
#' #' x = init(cna = x$clonal, ...)
#' #'
#' #'
#' #' # work with subclonal calls (omitting mutations and other parameters here)
#' #' x = analyze_peaks_subclonal(cna = x$subclonal, ...)
#' #' }
#' parse_Battenberg = function(x)
#' {
#'   required = c(
#'     "chr",
#'     "from",
#'     "to",
#'     "battenberg_nMaj1_A",
#'     "battenberg_nMin1_A",
#'     "battenberg_nMaj2_A",
#'     "battenberg_nMin2_A",
#'     "battenberg_frac1_A",
#'     "battenberg_frac2_A"
#'   )
#'
#'   if (!all(required %in% colnames(x)))
#'     stop(
#'       cli::format_error(
#'         "Missing Battenberg columns in these calls! Required columns are: {.field {required}}"
#'       )
#'     )
#'
#'   all_calls = x %>%
#'     dplyr::select(chr, from, to, starts_with("Battenberg")) %>%
#'     dplyr::rename(
#'       Major_1 = battenberg_nMaj1_A,
#'       minor_1 = battenberg_nMin1_A,
#'       Major_2 = battenberg_nMaj2_A,
#'       minor_2 = battenberg_nMin2_A,
#'       CCF_1 = battenberg_frac1_A,
#'       CCF_2 = battenberg_frac2_A
#'     ) %>%
#'     dplyr::select(chr, from, to, Major_1, minor_1, Major_2, minor_2, CCF_1, CCF_2) %>%
#'     dplyr::filter(!is.na(CCF_1))
#'
#'   clonal_calls = all_calls %>%
#'     dplyr::filter(CCF_1 == 1) %>%
#'     dplyr::rename(Major = Major_1,
#'                   minor = minor_1,
#'                   CCF = CCF_1) %>%
#'     dplyr::select(-Major_2,-minor_2,-CCF_2)
#'
#'   subclonal_calls = all_calls %>%
#'     dplyr::filter(!is.na(CCF_2))
#'
#'   return(list(clonal = clonal_calls, subclonal = subclonal_calls))
#' }
#'
#' x = readRDS("~/Downloads/8888e808-594b-4c76-b2e4-62aa56736f7c.rds")
#' x$cna %>% parse_Battenberg()
#'
#' x = init(x$mutations,)
#'
#
# x = readRDS("~/Downloads/f48c3c82-bebe-4b8e-909e-e1a51a7142ec.rds")
# x = init(x$mutations, x$cna, x$metadata$purity)
# x = x %>% analyze_peaks(n_bootstrap = 10)
# plot_peaks_analysis(x, what = 'general')
