##########################################
# Peak detection for simple karyotypes (more sofisticated algorithm) - rightmost peak matching
peak_detector = function(snvs,
                         expectation,
                         tumour_purity,
                         filtered_qc_snvs,
                         VAF_tolerance = 0.001,
                         matching_epsilon = 0.015,
                         kernel_adjust = 1,
                         p = 0.005,
                         KDE = FALSE)
{
  # expectation = CNAqc:::expected_vaf_peak(AB[1], AB[2], tumour_purity)
  # snvs = w %>%
  #   dplyr::sample_n(size = nrow(w), replace = TRUE)

  xy_peaks = den = NULL

  # Smoothed Gaussian kernel for VAF
  y = snvs %>% dplyr::pull(VAF)

  den = density(y, kernel = 'gaussian', adjust = kernel_adjust)
  in_range = den$x >= min(y) & den$x <= max(y)
  in_range = TRUE


  input_peakdetection = matrix(cbind(x = den$x[in_range], y = den$y[in_range]), ncol = 2)
  colnames(input_peakdetection) = c('x', 'y')

  # Test 5 parametrisations of peakPick neighlim
  pks = Reduce(dplyr::bind_rows,
               lapply(1:5,
                      function(n) {
                        pk = peakPick::peakpick(mat = input_peakdetection, neighlim = n)
                        input_peakdetection[pk[, 2], , drop = FALSE] %>% as.data.frame()
                      })) %>%
    as_tibble() %>%
    dplyr::arrange(x) %>%
    dplyr::mutate(x = round(x, 2), y = round(y, 2)) %>%
    dplyr::distinct(x, .keep_all = TRUE)


  # print(pks)
  hst = hist(snvs$VAF, breaks = seq(0, 1, 0.01), plot = F)$counts
  pks$counts_per_bin = hst[round(pks$x * 100)]

  # xy_peaks = pks %>%
  #   # dplyr::mutate(discarded = counts_per_bin < sum(hst) * p)
  #   dplyr::mutate(discarded = y <= 0.01)

  # Use KDE if required only
  if (KDE)
  {
    # Heuristic to remove low-density peaks
    xy_peaks = pks %>%
      # dplyr::mutate(discarded = counts_per_bin < sum(hst) * p)
      dplyr::mutate(discarded = y <= max(pks$y) * (1 / 20),
                    from = 'KDE')

    if (any(xy_peaks$discarded))
    {
      # cli::cli_alert_warning("Some KDE peaks have been removed")
      # xy_peaks %>%
      #   filter(discarded) %>%
      #   print
    }

  }

  # BMix clustering
  if (y %>% unique() %>% length() > 1)
  {
    bm = bmixfit(
      data.frame(successes = snvs$NV, trials = snvs$DP),
      K.BetaBinomials = 0,
      K.Binomials = 1:4,
      silent = TRUE
    )
    # plot.bmix(bm, data = data.frame(successes = rc$NV, trials = rc$DP))


    llxy = NULL
    for (b in names(bm$B.params))
    {
      w_den = which.min(abs(den$x - bm$B.params[b]))

      tnw = tibble(x = den$x[w_den],
                   y = den$y[w_den],
                   # counts_per_bin = bm$pi[b] * (snvs %>% nrow), # Wrong
                   discarded = FALSE)

      # Counts are counted the same way regardless it is a BMix fit or not.
      tnw$counts_per_bin = hst[round(tnw$x * 100)]

      llxy = llxy %>%
        bind_rows(tnw)
    }

    xy_peaks = xy_peaks %>% bind_rows(llxy %>% mutate(from = 'BMix'))
  }


  # Handle special case where everything is discarded by including the one
  # with highest value of counts_per_bin (just that).
  if (all(xy_peaks$discarded)) {
    xy_peaks = xy_peaks %>%
      mutate(discarded = ifelse(
        counts_per_bin == max(xy_peaks$counts_per_bin),
        TRUE,
        discarded
      ))
  }

  # Again, if the if-clause above did not suffice to find peaks, raise a stop error
  if (all(xy_peaks$discarded))
    stop("Cannot find peaks for this karyotype, raising an error (check the data distribution).")

  # linear combination of the weight, split by number of peaks to match
  weight = filtered_qc_snvs %>%
    filter(karyotype == snvs$karyotype[1]) %>%
    pull(norm_prop)

  weight = weight / nrow(expectation)

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
  if (nrow(match_xy_peaks) < nrow(expectation)) {
    entry = match_xy_peaks[1,]
    missing = nrow(expectation) - nrow(match_xy_peaks)
    for (s in 1:missing)
      match_xy_peaks = bind_rows(match_xy_peaks, entry)
  }

  expectation = dplyr::bind_cols(expectation, match_xy_peaks)

  # Distance in VAF space, converted to purity space
  expectation = expectation %>%
    mutate(
      offset_VAF = peak - x,
      # VAF space
      offset = compute_delta_purity(
        vaf = x,
        delta_vaf = offset_VAF,
        ploidy = strsplit(snvs$karyotype[1], ':') %>% unlist %>% as.numeric() %>% sum(),
        multiplicity = mutation_multiplicity
      )
    )

  # expectation$matched = abs(expectation$offset) <= matching_epsilon
  expectation$weight = weight
  expectation$epsilon = matching_epsilon
  expectation$VAF_tolerance = VAF_tolerance
  expectation$matched = NA

  for (i in  1:nrow(expectation))
    expectation$matched[i] = overlap_bands(
      peak = expectation$x[i],
      tolerance = VAF_tolerance,
      left_extremum = expectation$peak[i] - matching_epsilon[i],
      right_extremum = expectation$peak[i] + matching_epsilon[i]
    )


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

  return(list(
    matching = expectation,
    # plot = plot_data,
    density = den,
    xy_peaks = xy_peaks
  ))
}

##########################################
# Peak detection for simple karyotypes (more sofisticated algorithm) - closest peak matching
peak_detector_closest_hit_match = function(snvs,
                                           expectation,
                                           tumour_purity,
                                           filtered_qc_snvs,
                                           VAF_tolerance = 0.001,
                                           matching_epsilon = 0.015,
                                           kernel_adjust = 1,
                                           p = 0.005,
                                           KDE = FALSE)
{
  xy_peaks = den = NULL

  # Smoothed Gaussian kernel for VAF
  y = snvs %>% dplyr::pull(VAF)

  den = density(y,
                kernel = 'gaussian',
                adjust = kernel_adjust,
                na.rm = T)
  in_range = den$x >= min(y, na.rm = T) & den$x <= max(y, na.rm = T)

  # den = density(y, kernel = 'gaussian', adjust = 0.5)
  # plot(den)

  input_peakdetection = matrix(cbind(x = den$x[in_range], y = den$y[in_range]), ncol = 2)
  colnames(input_peakdetection) = c('x', 'y')

  # Test 5 parametrisations of peakPick neighlim
  pks = Reduce(dplyr::bind_rows,
               lapply(1:5,
                      function(n) {
                        pk = peakPick::peakpick(mat = input_peakdetection, neighlim = n)
                        input_peakdetection[pk[, 2], , drop = FALSE] %>% as.data.frame()
                      })) %>%
    as_tibble() %>%
    dplyr::arrange(x) %>%
    dplyr::mutate(x = round(x, 2), y = round(y, 2)) %>%
    dplyr::distinct(x, .keep_all = TRUE)

  # print(pks)
  hst = hist(snvs$VAF, breaks = seq(0, 1, 0.01), plot = F)$counts
  pks$counts_per_bin = hst[round(pks$x * 100)]

  # Run KDE if required only
  if (KDE)
  {
    # Heuristic to remove low-density peaks
    xy_peaks = pks %>%
      # dplyr::mutate(discarded = counts_per_bin < sum(hst) * p)
      dplyr::mutate(discarded = y <= max(pks$y) * (1 / 20),
                    from = 'KDE')


    if (any(xy_peaks$discarded))
    {
      # cli::cli_alert_warning("Some peaks have been removed")
      # xy_peaks %>%
      #   filter(discarded) %>%
      #   print
    }
  }

  # BMix clustering
  if (y %>% unique() %>% length() > 1)
  {
    # invisible(capture.output(
    bm = BMix::bmixfit(
      data.frame(successes = snvs$NV, trials = snvs$DP),
      K.BetaBinomials = 0,
      K.Binomials = 1:4,
      silent = TRUE
    )
    # ))
    # plot.bmix(bm, data = data.frame(successes = rc$NV, trials = rc$DP))

    llxy = NULL
    for (b in names(bm$B.params))
    {
      w_den = which.min(abs(den$x - bm$B.params[b]))

      tnw = tibble(x = den$x[w_den],
                   y = den$y[w_den],
                   # counts_per_bin = bm$pi[b] * (snvs %>% nrow), # Wrong
                   discarded = FALSE)

      # Counts are counted the same way regardless it is a BMix fit or not.
      tnw$counts_per_bin = hst[round(tnw$x * 100)]

      llxy = llxy %>%
        bind_rows(tnw)
    }

    xy_peaks = xy_peaks %>% bind_rows(llxy %>% mutate(from = 'BMix'))
  }

  # tibble(
  # `x`= bm$B.params,
  # `y` =
  # )


  # Handle special case where everything is discarded by including the one
  # with highest value of counts_per_bin (just that).
  if (all(xy_peaks$discarded)) {
    xy_peaks = xy_peaks %>%
      mutate(discarded = ifelse(
        counts_per_bin == max(xy_peaks$counts_per_bin),
        TRUE,
        discarded
      ))
  }

  # Again, if the if-clause above did not suffice to find peaks, raise a stop error
  if (all(xy_peaks$discarded))
    stop("Cannot find peaks for this karyotype, raising an error (check the data distribution).")


  # linear combination of the weight, split by number of peaks to match
  weight = filtered_qc_snvs %>%
    filter(karyotype == snvs$karyotype[1]) %>%
    pull(norm_prop)

  weight = weight / nrow(expectation)

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

    if (length(id_match) == 0)
      return(NA)
    else
      return(not_discarded[id_match,])
  }

  matched_peaks = lapply(seq_along(expectation$peak), get_match)

  # Matching table
  matching = expectation %>%
    dplyr::bind_cols(Reduce(dplyr::bind_rows, matched_peaks))

  # Distance in VAF space, converted to purity space
  matching = matching %>%
    rowwise() %>%
    mutate(
      offset_VAF = peak - x,
      # VAF space
      offset = compute_delta_purity(
        vaf = x,
        delta_vaf = offset_VAF,
        ploidy = strsplit(snvs$karyotype[1], ':') %>% unlist %>% as.numeric() %>% sum(),
        multiplicity = mutation_multiplicity
      )
    )

  # matching$matched = abs(matching$offset) <= matching_epsilon
  matching$weight = weight
  matching$epsilon = matching_epsilon

  matching$VAF_tolerance = VAF_tolerance
  matching$matched = NA

  for (i in  1:nrow(expectation))
    matching$matched[i] = overlap_bands(
      peak = matching$x[i],
      tolerance = VAF_tolerance,
      left_extremum = matching$peak[i] - matching_epsilon[i],
      right_extremum = matching$peak[i] + matching_epsilon[i]
    )

  # Density estimated
  density = den

  # Peaks
  peaks = xy_peaks

  return(list(
    matching = matching,
    density = den,
    xy_peaks = xy_peaks
  ))
}

##########################################
# Auxiliary functions  peak-matching for simple karyotypes
overlap_bands = function(peak,
                         tolerance,
                         left_extremum,
                         right_extremum)
{
  # if(peak + tolerance >= left_extremum & peak + tolerance <= right_extremum)
  #   return(TRUE)
  #
  # if(peak - tolerance >= left_extremum & peak - tolerance <= right_extremum)
  #   return(TRUE)

  M = max(peak - tolerance, left_extremum) <= min(peak + tolerance, right_extremum)

  return(M)
}

##########################################
# KDE-based pure peak detection (for general karyptypes and subclonal CNAs)
simple_peak_detector = function(mutations, kernel_adjust) {
  xy_peaks = den = NULL

  # Smoothed Gaussian kernel for VAF
  y = mutations %>% dplyr::pull(VAF)

  den = density(y,
                kernel = 'gaussian',
                adjust = kernel_adjust,
                na.rm = T)
  # in_range = den$x >= min(y, na.rm = T) & den$x <= max(y, na.rm = T)
  in_range = TRUE

  input_peakdetection = matrix(cbind(x = den$x[in_range], y = den$y[in_range]), ncol = 2)
  colnames(input_peakdetection) = c('x', 'y')

  # Test 5 parametrisations of peakPick neighlim
  pks = Reduce(dplyr::bind_rows,
               lapply(1:5,
                      function(n) {
                        pk = peakPick::peakpick(mat = input_peakdetection, neighlim = n)
                        input_peakdetection[pk[, 2], , drop = FALSE] %>% as.data.frame()
                      })) %>%
    as_tibble() %>%
    dplyr::arrange(x) %>%
    dplyr::mutate(x = round(x, 2), y = round(y, 2)) %>%
    dplyr::distinct(x, .keep_all = TRUE) %>%
    dplyr::mutate(x = case_when(x > 1 & x < 1.01 ~ 1,
                                x < 0 & x > -0.01 ~ 0,
                                TRUE ~ x),) %>%
    dplyr::filter(x <= 1, x >= 0)

  hst = hist(mutations$VAF,
             breaks = seq(0, 1, 0.01),
             plot = F)$counts
  pks$counts_per_bin = hst[round(pks$x * 100)]

  # Heuristic to remove low-density peaks
  pks = pks %>%
    dplyr::mutate(discarded = y <= max(pks$y) * (1 / 20), from = 'KDE')

  return(list(peaks = pks, density = den))
}
