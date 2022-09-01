



# data('example_dataset_CNAqc', package = 'CNAqc')
# x = init(example_dataset_CNAqc$mutations, example_dataset_CNAqc$cna,example_dataset_CNAqc$purity)



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
