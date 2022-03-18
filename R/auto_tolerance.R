#' Determine the optimal error tolerance based on data
#'
#' @description Regression has been used to measure the rate of false positives
#' from simulated tumours with variable coverage and purity. This allows to determine
#' an optimal valut of $\epsilon$, parameter `purity_error` of function
#' \code{analyze_peaks}.
#'
#' @param purity Data purity (putative).
#' @param coverage Data coverage.
#' @param fpr Desired false positive rate.
#' @param epsilon_range Range of values to constrain $epsilon$.
#'
#' @return The $\epsilon$ value estimated from data, constrained to be in
#' `epsilon_range`, in order to limit the false positive rate to be `fpr`.
#'
#' @importFrom  gtools mixedsort
#'
#' @export
#'
#' @examples
#' auto_tolerance(.3, 90)
auto_tolerance = function(
  purity,
  coverage,
  fpr = 0.1,
  epsilon_range = c(0.01, 0.08)
)
{
  fpr_test = CNAqc::fpr_test %>% filter(purity != .60)

  fake_ret = list(
    epsilon_range = 0.01,
    training_plot = ggplot(),
    extrapolation_plot = ggplot()
  )

  # Outside training set
  if(fpr_test$coverage %>% max < coverage)
  {
    cli::cli_alert_danger("Coverage {.field {coverage}} is unsupported - maximum is {.value {fpr_test$coverage %>% max}} - returning 0.01 (default)")
    return(fake_ret)
  }

  if(fpr_test$coverage %>% min > coverage)
  {
    cli::cli_alert_danger("Coverage {.field {coverage}} is unsupported - minimum is {.value {fpr_test$coverage %>% min}} - returning 0.01 (default)")
    return(fake_ret)
  }

  if(fpr_test$purity %>% max < purity)
  {
    cli::cli_alert_danger("Purity {.field {purity}} is unsupported - maximum is {.value {fpr_test$purity %>% max}} - returning 0.01 (default)")
    return(fake_ret)
  }

  if(fpr_test$purity %>% min > purity)
  {
    cli::cli_alert_danger("Purity {.field {purity}} is unsupported - minimum is {.value {fpr_test$purity %>% min}} - returning 0.01 (default)")
    return(fake_ret)
  }

  if(fpr >= 1 | fpr <= 0)
  {
    cli::cli_alert_danger("FPR outside [0,1] does not make sense - returning 0.01 (default)")
    return(fake_ret)
  }

  # We can work and regress the input
  fpr_test = fpr_test %>%
    dplyr::mutate(key = paste0(coverage, ":", purity))

  # fpr_test %>%
  #   filter(coverage == 60) %>%
  #   ggplot(aes(x=epsilon_tolerance, y = FPR, color = purity %>% paste())) +
  #   geom_point() +
  #   geom_smooth(method = 'lm')

  # Stat smooth data - create linear functions
  ggp_data = fpr_test %>%
    dplyr::group_split(key) %>%
    lapply(function(w)
    {
      fit  <-
        glm(FPR ~ epsilon_tolerance,
            data = w %>%  dplyr::select(epsilon_tolerance, FPR))

      # y = mx + q
      q = fit$coefficients[1] %>% as.numeric()
      m = fit$coefficients[2] %>% as.numeric()

      # plot(w$epsilon_tolerance, w$FPR)
      # points(w$epsilon_tolerance, m * (w$epsilon_tolerance) + q, col = 'red', pch = 11)

      # x = (y - q)/m
      inv_fun = function(FPR) {
        eval = (FPR - q) / m
        if (eval < epsilon_range[1]) {
          warning(
            "Tolerance below 1% for FPR = ",
            FPR,
            " - returning ",
            epsilon_range[1],
            " (epsilon_range)."
          )
          return(epsilon_range[1])
        }

        if (eval > epsilon_range[2]) {
          warning(
            "Tolerance above 100% for FPR = ",
            FPR,
            " - returning ",
            epsilon_range[1],
            " (epsilon_range)."
          )
          return(epsilon_range[2])
        }

        eval
      }

      w[1,] %>%
        dplyr::select(coverage, purity) %>%
        dplyr::mutate(epsilon = inv_fun(fpr))
    })

  regression_data = Reduce(f = bind_rows, ggp_data)

  # names(ggp_data) = fpr_test$key %>% unique

  # Evaluate regression test - inverse function
  cli::cli_alert("Inverting training from regression; requireing epsilon in range {.field [{epsilon_range[1]} - {epsilon_range[2]}]}.")

  # regression_data = tidyr::expand_grid(coverage = fpr_test$coverage %>% unique,
  #                                      purity = fpr_test$purity %>% unique) %>%
  #   mutate(fun = ggp_data[paste0(coverage, ':', purity)])
  #
  # regression_data$fun_x = sapply(regression_data$fun, function(f)
  #   f(fpr) %>% round(4))

  # Tile plot
  reg_plot2 = regression_data %>%
    dplyr::mutate(coverage = factor(
      coverage,
      levels = gtools::mixedsort(regression_data$coverage %>% unique)
    )) %>%
    mutate(purity = factor(
      purity,
      levels = gtools::mixedsort(regression_data$purity %>% unique)
    )) %>%
    ggplot() +
    geom_tile(aes(
      x = coverage,
      y = purity,
      fill = epsilon,
      width = 0.9,
      height = 0.9
    )) +
    CNAqc:::my_ggplot_theme() +
    scale_fill_gradientn(colours = c('indianred3', 'forestgreen', 'steelblue')) +
    # scale_fill_viridis_c(option = 'C', limits = c(0, NA)) +
    guides(fill = guide_colorbar(bquote("Regressed " * epsilon * "  "), barwidth = unit(3, 'cm'))) +
    labs(title = bquote("Regressed " * epsilon * ' for FPR < ' * .(fpr)))

  # 2D interpolation
  cli::cli_alert("Generating interpolated 2D regression map.")

  purity_x = seq(min(fpr_test$purity), max(fpr_test$purity), 0.01)
  coverage_x = seq(min(fpr_test$coverage), max(fpr_test$coverage), 1)

  plot_grid = tidyr::expand_grid(coverage = coverage_x,
                                 purity = purity_x) %>%
    rowwise() %>%
    mutate(
      fit =
        akima::interp(
          regression_data$purity,
          regression_data$coverage,
          regression_data$epsilon,
          xo = purity,
          yo = coverage,
          linear = TRUE,
          extrap = TRUE
        )$z %>% as.numeric
    )

  # Shape data, apply cutoff
  # cli::cli_alert("Generating interpolated regression map.")

  plot_grid = plot_grid %>%
    filter(!is.na(fit)) %>%
    mutate(coverage = factor(coverage,
                             levels = gtools::mixedsort(plot_grid$coverage %>% unique))) %>%
    mutate(purity = factor(purity,
                           levels = gtools::mixedsort(plot_grid$purity %>% unique))) %>%
    mutate(fit = ifelse(fit > epsilon_range[2], epsilon_range[2], fit)) %>%
    mutate(fit = ifelse(fit < epsilon_range[1], epsilon_range[1], fit)) %>%
    ggplot() +
    geom_tile(aes(x = coverage, y = purity, fill = fit)) +
    scale_y_discrete(breaks = regression_data$purity %>% unique) +
    CNAqc:::my_ggplot_theme() +
    # scale_fill_viridis_c(option = 'C', limits = c(0, NA)) +
    guides(fill = guide_colorbar(bquote("Extrapolated " * .(epsilon_range[1]) * "<" * epsilon * "<" * .(epsilon_range[2]) *"  "), barwidth = unit(3, 'cm'))) +
    labs(title = bquote("Extrapolated " * epsilon * ' map')) +
    scale_fill_gradientn(colours = c('indianred3', 'forestgreen', 'steelblue')) +
    scale_x_discrete(breaks = regression_data$coverage %>% unique)

  # Point_estimate
  point_estimate = akima::interp(
    regression_data$purity,
    regression_data$coverage,
    regression_data$epsilon,
    xo = purity,
    yo = coverage,
    linear = TRUE,
    extrap = TRUE
  )$z %>% as.numeric

  if(point_estimate < epsilon_range[1]) point_estimate = epsilon_range[1]
  if(point_estimate > epsilon_range[2]) point_estimate = epsilon_range[2]

  p_est = data.frame(coverage = coverage, purity = purity) %>%
    mutate(coverage = factor(coverage,
                             levels = gtools::mixedsort(coverage_x))) %>%
    mutate(purity = factor(purity,
                           levels = gtools::mixedsort(purity_x)))

  plot_grid = plot_grid +
    geom_point(
      data = p_est,
      aes(x = coverage, y = purity),
      fill = 'black',
      size = 2,
      pch = 23
    )

  cli::cli_alert_success("Suggested point estimate: {.field \u03B5 = {point_estimate}}.")

  return(
    list(
      epsilon_range = point_estimate,
      training_plot = reg_plot2,
      extrapolation_plot = plot_grid
    )
  )
}

I
