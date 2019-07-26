
peak_detector = function(snvs,
                         expectation,
                         tumour_purity,
                         filtered_qc_snvs,
                         adjust = 1,
                         neighlim = 1,
                         min_density = 0.5,
                         matching_epsilon = 0.015,
                         ...)
{
  y = snvs %>% pull(VAF)
  den = density(y, kernel = 'gaussian', adjust = adjust)

  input_peakdetection = matrix(cbind(x = den$x, y = den$y), ncol = 2)
  colnames(input_peakdetection) = c('x', 'y')

 # pio::pioStr("peak_detector", 'karyotype', snvs$karyotype[1], ' - peakPick params ...')
  require(peakPick)
  peaks =  peakpick(mat = input_peakdetection, neighlim = neighlim, ...)
  peaks_params = list(...)

  xy_peaks = input_peakdetection[peaks[, 2], , drop = FALSE] %>%
    as_tibble() %>%
    mutate(
      x = round(x, 2),
      y = round(y, 2)
    ) %>%
    mutate(discarded = y < min_density)

  # linear combination of the weight, split by number of peaks to match
  weight = filtered_qc_snvs %>%
    filter(karyotype == snvs$karyotype[1]) %>%
    pull(norm_prop)

  weight = weight/nrow(expectation)


  # Plot
  plot_data = ggplot(data = snvs, aes(VAF)) +
    geom_histogram(aes(y = ..density..), binwidth = 0.01, alpha = .3) +
    geom_line(
      data = data.frame(x = den$x, y = den$y),
      aes(x = x, y = y),
      size = .3,
      color = 'black'
    ) +
    theme_minimal() +
    labs(
      title = paste0("Karyotype ", snvs$karyotype[1]),
      subtitle = paste0('w = ', round(weight, 3), ' (n = ', nrow(snvs), ')'),
      # caption = paste0(
      #   "  KDE bandwith adjustment : ", adjust, '\n',
      #   "peakPick neighlim integer : ", neighlim
      # ),
      y = 'KDE'
    ) +
    geom_hline(
      yintercept = min_density,
      color = 'darkred',
      linetype = 'dashed',
      size = .3
    ) +
    theme(
      legend.position = 'bottom'
    )  +
    xlim(0, 1)

  # Add points for peaks to plot
  plot_data = plot_data +
    geom_point(data = xy_peaks,
               aes(x = x, y = y, shape = discarded),
               size = 1.5,
               show.legend = FALSE
    ) +
    scale_shape_manual(values = c(`TRUE` = 1, `FALSE` = 16))

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

  expectation = bind_cols(expectation, match_xy_peaks) %>%
    mutate(
      offset = peak - x,
      matched = abs(offset) <= matching_epsilon
    )

  expectation$weight = weight

  # Add expectation peaks, and matching colors
  plot_data = plot_data +
    geom_point(
      data = expectation,
      aes(x = x, y = y, color = matched),
      size = 2,
      show.legend = FALSE
    ) +
    geom_segment(
      data = expectation,
      aes(x = x, y = y, xend = peak, yend = y, color = matched),
      show.legend = FALSE
    )+
    annotate(
      geom = 'rect',
      xmin = expectation$peak - matching_epsilon,
      xmax = expectation$peak + matching_epsilon,
      ymin = 0,
      ymax = Inf,
      color = NA,
      alpha = .4,
      fill = 'steelblue'
    ) +
    geom_vline(
      data = expectation,
      aes(xintercept = peak, color = matched),
      size = .7,
      linetype = 'longdash',
      show.legend = FALSE
    ) +
    scale_color_manual(values = c(`TRUE` = 'forestgreen', `FALSE` = 'red'))

  # Annotate the offset number
  plot_data = plot_data +
  ggrepel::geom_text_repel(
    data = expectation %>% filter(!matched),
    aes(x = x, y = y, label = round(offset, 2), color = matched),
    size = 3,
    nudge_x = 0,
    nudge_y = 0,
    show.legend = FALSE
  )

  return(list(matching = expectation, plot = plot_data))
}



