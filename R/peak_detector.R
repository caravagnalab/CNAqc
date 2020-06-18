
peak_detector = function(snvs,
                         expectation,
                         tumour_purity,
                         filtered_qc_snvs,
                         matching_epsilon = 0.015,
                         kernel_adjust = 1,
                         p = 0.005)
{
  # expectation = CNAqc:::expected_vaf_peak(AB[1], AB[2], tumour_purity)
  # snvs = w %>%
  #   dplyr::sample_n(size = nrow(w), replace = TRUE)

  # Smoothed Gaussian kernel for VAF
  y = snvs %>% dplyr::pull(VAF)

  den = density(y, kernel = 'gaussian', adjust = kernel_adjust)
  in_range = den$x >= min(y) & den$x <= max(y)

  input_peakdetection = matrix(cbind(x = den$x[in_range], y = den$y[in_range]), ncol = 2)
  colnames(input_peakdetection) = c('x', 'y')

  # Test 5 parametrisations of peakPick neighlim
  pks = Reduce(
    dplyr::bind_rows,
    lapply(1:5,
           function(n){
             pk = peakPick::peakpick(mat = input_peakdetection, neighlim = n)
             input_peakdetection[pk[, 2], , drop = FALSE] %>% as.data.frame()
           }
    )
  ) %>%
    as_tibble() %>%
    dplyr::arrange(x) %>%
    dplyr::mutate(x = round(x, 2), y = round(y, 2)) %>%
    dplyr::distinct(x, .keep_all = TRUE)


  # print(pks)
  hst = hist(snvs$VAF, breaks = seq(0, 1, 0.01), plot = F)$counts
  pks$counts_per_bin = hst[round(pks$x * 100)]

  xy_peaks = pks %>%
    dplyr::mutate(discarded = counts_per_bin < sum(hst) * p)

  # linear combination of the weight, split by number of peaks to match
  weight = filtered_qc_snvs %>%
    filter(karyotype == snvs$karyotype[1]) %>%
    pull(norm_prop)

  weight = weight/nrow(expectation)

  ###### ###### ###### ###### ######
  # # Plot
  # plot_data = ggplot(data = snvs, aes(VAF)) +
  #   geom_histogram(aes(y = ..density..), binwidth = 0.01, alpha = .3) +
  #   geom_line(
  #     data = data.frame(x = den$x, y = den$y),
  #     aes(x = x, y = y),
  #     size = .3,
  #     color = 'black'
  #   ) +
  #   CNAqc:::my_ggplot_theme() +
  #   labs(
  #     title = paste0("Karyotype ", snvs$karyotype[1]),
  #     subtitle = paste0('w = ', round(weight, 3), ' (n = ', nrow(snvs), ')'),
  #     y = 'KDE'
  #   ) +
  #   # geom_hline(
  #   #   yintercept = sum(hst) * 0.005,
  #   #   color = 'darkred',
  #   #   linetype = 'dashed',
  #   #   size = .3
  #   # ) +
  #   theme(
  #     legend.position = 'bottom'
  #   )  +
  #   xlim(0, 1)

  # Add points for peaks to plot
  # plot_data = plot_data +
  #   geom_point(data = xy_peaks,
  #              aes(x = x, y = y, shape = discarded, size = counts_per_bin),
  #              # size = 1.5,
  #              show.legend = FALSE
  #   ) +
  #   scale_shape_manual(values = c(`TRUE` = 1, `FALSE` = 16)) +
  #   scale_size(range = c(1, 3))

  # Compute matching
  #
  # 1) sort the expected peaks and the xy_peaks by expected VAF
  # 2) combine them and subtract
  expectation = expectation %>% arrange(desc(peak))
  match_xy_peaks = xy_peaks %>%
    filter(!discarded) %>%
    arrange(desc(x)) %>%
    filter(row_number() <= nrow(expectation))

  # Handle the special case where we have less peaks that the ones we need to
  # match. In this case we match everything to the same peak
  if(nrow(match_xy_peaks) < nrow(expectation)) {
    entry = match_xy_peaks[1, ]
    missing = nrow(expectation) - nrow(match_xy_peaks)
    for(s in 1:missing) match_xy_peaks = bind_rows(match_xy_peaks, entry)
  }

  expectation = dplyr::bind_cols(expectation, match_xy_peaks) %>%
    dplyr::mutate(
      offset = peak - x,
      matched = abs(offset) <= matching_epsilon
    )

  expectation$weight = weight

  # Add expectation peaks, and matching colors
  # plot_data = plot_data +
  #   geom_point(
  #     data = expectation,
  #     aes(x = x, y = y, color = matched),
  #     size = 2,
  #     shape = 4,
  #     show.legend = FALSE
  #   ) +
  #   geom_segment(
  #     data = expectation,
  #     aes(x = x, y = y, xend = peak, yend = y, color = matched),
  #     show.legend = FALSE
  #   )+
  #   annotate(
  #     geom = 'rect',
  #     xmin = expectation$peak - matching_epsilon,
  #     xmax = expectation$peak + matching_epsilon,
  #     ymin = 0,
  #     ymax = Inf,
  #     color = NA,
  #     alpha = .4,
  #     fill = 'steelblue'
  #   ) +
  #   geom_vline(
  #     data = expectation,
  #     aes(xintercept = peak, color = matched),
  #     size = .7,
  #     linetype = 'longdash',
  #     show.legend = FALSE
  #   ) +
  #   scale_color_manual(values = c(`TRUE` = 'forestgreen', `FALSE` = 'red'))

  # Annotate the offset number
  # cex_opt = getOption('CNAqc_cex', default = 1)

  # plot_data = plot_data +
  # ggrepel::geom_text_repel(
  #   data = expectation %>% filter(!matched),
  #   aes(x = x, y = y, label = round(offset, 2), color = matched),
  #   nudge_x = 0,
  #   nudge_y = 0,
  #   size = 3 * cex_opt,
  #   show.legend = FALSE
  # )
  # # plot_data

  return(
    list(
    matching = expectation,
    # plot = plot_data,
    density = den,
    xy_peaks = xy_peaks)
    )
}

peak_detector_closest_hit_match = function(snvs,
                         expectation,
                         tumour_purity,
                         filtered_qc_snvs,
                         matching_epsilon = 0.015,
                         kernel_adjust = 1,
                         p = 0.005)
{
  # Smoothed Gaussian kernel for VAF
  y = snvs %>% dplyr::pull(VAF)

  den = density(y, kernel = 'gaussian', adjust = kernel_adjust)
  in_range = den$x >= min(y) & den$x <= max(y)

  input_peakdetection = matrix(cbind(x = den$x[in_range], y = den$y[in_range]), ncol = 2)
  colnames(input_peakdetection) = c('x', 'y')

  # Test 5 parametrisations of peakPick neighlim
  pks = Reduce(
    dplyr::bind_rows,
    lapply(1:5,
           function(n){
             pk = peakPick::peakpick(mat = input_peakdetection, neighlim = n)
             input_peakdetection[pk[, 2], , drop = FALSE] %>% as.data.frame()
           }
    )
  ) %>%
    as_tibble() %>%
    dplyr::arrange(x) %>%
    dplyr::mutate(x = round(x, 2), y = round(y, 2)) %>%
    dplyr::distinct(x, .keep_all = TRUE)


  # print(pks)
  hst = hist(snvs$VAF, breaks = seq(0, 1, 0.01), plot = F)$counts
  pks$counts_per_bin = hst[round(pks$x * 100)]

  xy_peaks = pks %>%
    dplyr::mutate(discarded = counts_per_bin < sum(hst) * p)

  # linear combination of the weight, split by number of peaks to match
  weight = filtered_qc_snvs %>%
    filter(karyotype == snvs$karyotype[1]) %>%
    pull(norm_prop)

  weight = weight/nrow(expectation)

  ###### ###### ###### ###### ######
  # Compute matching ~ get any possible match given the expectation
  #
  # 1) sort the expected peaks and the xy_peaks by expected VAF
  # 2) combine them and subtract
  get_match = function(p)
  {
    not_discarded = xy_peaks %>% filter(!discarded)

    distances = abs(not_discarded$x - expectation$peak[p])
    id_match = which.min(distances)

    if(length(id_match) == 0) return(NA)
    else return(not_discarded[id_match, ])
  }

  matched_peaks = lapply(seq_along(expectation$peak), get_match)

  # Matching table
  matching = expectation %>%
    dplyr::bind_cols(Reduce(dplyr::bind_rows, matched_peaks))

  matching$offset = matching$peak - matching$x
  matching$matched = matching$offset <= matching_epsilon
  matching$weight = weight

  # Density estimated
  density = den

  # Peaks
  peaks = xy_peaks

  return(
    list(
      matching = matching,
      density = den,
      xy_peaks = xy_peaks)
  )
}

