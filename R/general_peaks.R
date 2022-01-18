



# data('example_dataset_CNAqc', package = 'CNAqc')
# x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna,example_dataset_CNAqc$purity)

analyze_peaks_general = function(x,
                                 n_min = 50,
                                 epsilon = 0.03,
                                 kernel_adjust = 1,
                                 n_bootstrap = 5)
{
  # aux_peak_calling_general = function(mutations){
  #   xy_peaks = den = NULL
  #
  #   # Smoothed Gaussian kernel for VAF
  #   y = mutations %>% dplyr::pull(VAF)
  #
  #   den = density(y, kernel = 'gaussian', adjust = kernel_adjust, na.rm = T)
  #   in_range = den$x >= min(y, na.rm = T) & den$x <= max(y, na.rm = T)
  #
  #   input_peakdetection = matrix(cbind(x = den$x[in_range], y = den$y[in_range]), ncol = 2)
  #   colnames(input_peakdetection) = c('x', 'y')
  #
  #   # Test 5 parametrisations of peakPick neighlim
  #   pks = Reduce(dplyr::bind_rows,
  #                lapply(1:5,
  #                       function(n) {
  #                         pk = peakPick::peakpick(mat = input_peakdetection, neighlim = n)
  #                         input_peakdetection[pk[, 2], , drop = FALSE] %>% as.data.frame()
  #                       })) %>%
  #     as_tibble() %>%
  #     dplyr::arrange(x) %>%
  #     dplyr::mutate(x = round(x, 2), y = round(y, 2)) %>%
  #     dplyr::distinct(x, .keep_all = TRUE)
  #
  #   hst = hist(mutations$VAF, breaks = seq(0, 1, 0.01), plot = F)$counts
  #   pks$counts_per_bin = hst[round(pks$x * 100)]
  #
  #   # Heuristic to remove low-density peaks
  #   pks = pks %>%
  #     dplyr::mutate(discarded = y <= max(pks$y) * (1 / 20), from = 'KDE')
  #
  #   return(list(peaks = pks, density = den))
  # }

  # Filter small segments
  candidates = x$snvs$karyotype %>% unique
  candidates = setdiff(candidates,  c("1:1", "1:0", "2:0", "2:1", "2:2", "NA:NA"))
  n_cand = x$n_karyotype[candidates] >= n_min

  analysis = names(n_cand)[n_cand]

  # Expected peaks
  expected_peaks = lapply(analysis,
                          function(k)
                            expectations_generalised(
                              m = NULL,
                              M = NULL,
                              p = x$purity,
                              karyotype = k
                            )) %>%
    Reduce(f = bind_rows)

  # Data peaks and densities
  data_fit = lapply(analysis,
                    function(k) {
                      # Get default: all data run
                      s_run = simple_peak_detector(x$snvs %>% filter(karyotype == k), kernel_adjust = kernel_adjust)

                      if (n_bootstrap > 1)
                      {
                        # Bootstrap if required (and merge peaks)
                        N = x$snvs %>% filter(karyotype == k) %>% nrow()

                        sb_run = lapply(1:n_bootstrap, function(e) {
                          simple_peak_detector(
                            x$snvs %>%
                              filter(karyotype == k) %>%
                              sample_n(N, replace = TRUE),
                            kernel_adjust = kernel_adjust
                          )$peaks
                        }) %>%
                          Reduce(f = bind_rows) %>%
                          distinct(x, .keep_all = TRUE)

                        s_run$peaks = s_run$peaks %>%
                          bind_rows(sb_run) %>%
                          distinct(x, .keep_all = TRUE)
                      }
                      return(s_run)
                    })

  data_densities = lapply(data_fit %>% seq_along,
                          function(f) {
                            data.frame(
                              x = data_fit[[f]]$density$x,
                              y = data_fit[[f]]$density$y,
                              karyotype = analysis[f]
                            )
                          }) %>%
    Reduce(f = bind_rows)

  data_peaks = lapply(data_fit %>% seq_along,
                      function(f) {
                        data_fit[[f]]$peaks %>%
                          mutate(karyotype = analysis[f])
                      }) %>%
    Reduce(f = bind_rows)

  # Matching
  for (e in 1:nrow(expected_peaks))
  {
    e_k = expected_peaks$karyotype[e]
    d_k = abs(data_peaks %>%
                filter(karyotype == e_k) %>%
                pull(x) -
                expected_peaks$peak[e])

    if (any(d_k < epsilon))
      expected_peaks$matched[e] = TRUE
    else
      expected_peaks$matched[e] = FALSE
  }

  add_counts = function(w) {
    w$n = x$n_karyotype[w$karyotype]
    w %>% as_tibble()
  }

  # Results
  x$peaks_analysis$general$analysis = analysis
  x$peaks_analysis$general$params = list(n_min = n_min,
                                         epsilon = epsilon,
                                         n_bootstrap = n_bootstrap)
  x$peaks_analysis$general$expected_peaks = expected_peaks %>% add_counts
  x$peaks_analysis$general$data_peaks = data_peaks %>% add_counts
  x$peaks_analysis$general$data_densities = data_densities %>% add_counts

  # Summary table
  stable =
    x$peaks_analysis$general$expected_peaks %>%
    group_by(karyotype, matched, n) %>%
    summarise(hits = n()) %>%
    mutate(matched = ifelse(matched, 'matched', 'mismatched')) %>%
    pivot_wider(names_from = 'matched', values_from = 'hits')

  if ("matched" %in% (stable %>% colnames))
    stable = stable %>%
    mutate(matched = ifelse(is.na(matched), 0, matched))
  else
    stable$matched = 0

  if ("mismatched" %in% (stable %>% colnames))
    stable = stable %>%
    mutate(mismatched = ifelse(is.na(mismatched), 0, mismatched))
  else
    stable$mismatched = 0


  stable = stable %>%
    mutate(prop = matched / (matched + mismatched)) %>%
    arrange(prop %>% desc)

  x$peaks_analysis$general$summary = stable

  return(x)
}

# eps_seq = seq(0.01, .2, 0.01)
# crv = sapply(eps_seq, function(e){
#   w = analyze_peaks_general(x, epsilon = e, n_bootstrap = 1)
#   m = w$peaks_analysis$general$summary$matched %>% sum
#   mm = w$peaks_analysis$general$summary$mismatched %>% sum
#   m/(m+mm)
# })
#
# ggplot(data.frame(epsilon = eps_seq, proportion = crv)) +
#   geom_line(aes(x = epsilon, y = 1-proportion)) +
#   geom_point(aes(x = epsilon, y = 1-proportion))

auto_tune_epsilon = function(x,
                             epsilon = seq(0.01, .2, 0.01),
                             n_bootstrap = 1)
{
  cli::cli_h1("Autotune epsilon via peak matching")

  # Values for peaks
  values = easypar::run(
    FUN = function(e) {
      w = analyze_peaks_general(x, epsilon = e, n_bootstrap = n_bootstrap)

      if ('matched' %in% (w$peaks_analysis$general$summary %>% colnames))
        m = w$peaks_analysis$general$summary$matched %>% sum
      else
        m = 0

      if ('mismatched' %in% (w$peaks_analysis$general$summary %>% colnames))
        mm = w$peaks_analysis$general$summary$mismatched %>% sum
      else
        mm = 0

      return(m / (m + mm))
    },
    PARAMS = lapply(epsilon, list),
    parallel = FALSE
  )

  values = values %>% unlist()

  # Compute score
  score <- diff(values) / diff(epsilon) # first derivative
  score = abs(score) * values[-1]       # weighted by % passes
  score = score / max(score)              # rescaled to cap at 1

  score_max = which.max(score) + 2      # Best epsilon and best matches
  best_epsilon = epsilon[score_max]
  best_matched = round(values[score_max] * 100, 2)

  elbow_plot = ggplot(data.frame(epsilon = epsilon, proportion = values)) +
    geom_point(
      data = data.frame(x = best_epsilon, y = 1 - best_matched),
      aes(x = x, y = y),
      color = 'steelblue',
      fill = 'gray',
      alpha = 1,
      pch = 21,
      size = 5
    ) +
    geom_line(aes(x = epsilon, y = 1 - proportion)) +
    geom_point(aes(x = epsilon, y = 1 - proportion)) +
    CNAqc:::my_ggplot_theme() +
    geom_line(
      data = data.frame(x = epsilon[-1], y = d1),
      aes(x = x, y = y),
      color = 'red',
      linetype = 'dashed',
      size = .2
    ) +
    geom_point(
      data = data.frame(x = epsilon[-1], y = d1),
      aes(x = x, y = y),
      color = 'indianred3',
      size = 1,
      pch = 23
    ) +
    labs(
      y = "% of mismatches",
      x = bquote(epsilon),
      title = bquote("Auto-tune " * epsilon["opt"] * ' = ' * .(best_epsilon))
    )

  # Best
  cli::cli_alert_success(
    "Best epsilon estimated: \u03b5 = {.field {best_epsilon}} [% matched {.field {best_matched}%}]"
  )

  cli::cli_h2("Running peak-detection with best epsilon \u03b5 = {.field {best_epsilon}}")

  x = analyze_peaks(x, purity_error = best_epsilon)

  plot_qc_peaks = plot_peaks_analysis(x, what = 'common')
  plot_general_peaks = plot_peaks_analysis(x, what = 'general')

  cowplot::plot_grid(
    elbow_plot,
    pk_plot,
    rel_widths = c(1, 2),
    ncol = 2,
    align = 'h',
    axis = 'bt'
  )
}
