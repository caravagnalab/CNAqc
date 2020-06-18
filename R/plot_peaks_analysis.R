#' Plot the results of peak analysis.
#'
#'  @description Results from \code{analyze_peaks} can be visualised with this
#'  function, which arranges the plots of each karyotype in a figure via \code{ggpubr}.
#'  Each karyotype shows the data, the estimated density, the peaks (selected and
#'  discarded), and the fit with shaded matching area.
#'
#' @param x An object of class \code{cnaqc}, where function \code{analyze_peaks} has
#' been computed.
#'
#' @return A \code{ggpubr} object for an assembled figure.
#' @export
#'
#' @import ggpubr
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna,example_dataset_CNAqc$purity)
#'
#' x = analyze_peaks(x)
#' plot_peaks_analysis(x)
plot_peaks_analysis = function(x)
{
  stopifnot(inherits(x, "cnaqc"))

  with_peaks = all(!is.null(x$peaks_analysis))
  if (!with_peaks) {
    warning("Input does not have peaks, see ?peaks_analysis to run peaks analysis.")
    return(CNAqc:::eplot())
  }

  karyotypes = x$peaks_analysis$fits %>% names

  # Plot each one of the fits
  plots = lapply(karyotypes, function(k) {
    if (all(is.null(x$peaks_analysis$fits[[k]]$matching)))
      return(CNAqc:::eplot())

    return(suppressWarnings(suppressMessages(plot_peaks_fit(x, k))))
  })

  # Overall QC
  qc = ifelse(x$peaks_analysis$QC == 'PASS', 'forestgreen', 'indianred3')

  # Plots assembly
  suppressWarnings(suppressMessages(
    ggpubr::ggarrange(
      plotlist = plots,
      nrow = 1,
      ncol = length(plots)
    ) +
      theme(
        title = element_text(color = qc),
        panel.border = element_rect(colour = qc,
                                    fill = NA)
      )
  ))
}

# Plot a single run results
plot_peaks_fit = function(x, k)
{
  # Required input values
  snvs = x$snvs %>% dplyr::filter(karyotype == k)
  den = x$peaks_analysis$fits[[k]]$density
  expectation = x$peaks_analysis$fits[[k]]$matching
  xy_peaks = x$peaks_analysis$fits[[k]]$xy_peaks
  matching_epsilon = x$peaks_analysis$matching_epsilon

  # linear combination of the weight, split by number of peaks to match
  weight = x$peaks_analysis$matches %>%
    dplyr::filter(karyotype == k) %>%
    dplyr::pull(weight) %>%
    sum

  # Plots cex for anything that is not the main theme
  cex_opt = getOption('CNAqc_cex', default = 1)

  # Plot the data
  plot_data =
    ggplot(data = snvs, aes(VAF)) +
    geom_histogram(aes(y = ..density..), binwidth = 0.01, alpha = .3) +
    geom_line(
      data = data.frame(x = den$x, y = den$y),
      aes(x = x, y = y),
      size = .3,
      color = 'black'
    ) +
    CNAqc:::my_ggplot_theme() +
    labs(
      title = paste0("Karyotype ", snvs$karyotype[1]),
      subtitle = paste0('w = ', round(weight, 3), ' (n = ', nrow(snvs), ')'),
      y = 'KDE',
      x = paste0("VAF (match: ", x$peaks_analysis$matching_strategy, ")")
    ) +
    theme(legend.position = 'bottom')  +
    xlim(0, 1)

  # Add points for peaks to plot
  plot_data = plot_data +
    geom_point(
      data = xy_peaks,
      aes(
        x = x,
        y = y,
        shape = discarded,
        size = counts_per_bin
      ),
      show.legend = FALSE
    ) +
    scale_shape_manual(values = c(`TRUE` = 1, `FALSE` = 16)) +
    scale_size(range = c(1, 3) * cex_opt)

  # Add expectation peaks, and matching colors
  plot_data = plot_data +
    geom_point(
      data = expectation,
      aes(x = x, y = y, color = matched),
      size = 2 * cex_opt,
      shape = 4,
      show.legend = FALSE
    ) +
    geom_segment(
      data = expectation,
      aes(
        x = x,
        y = y,
        xend = peak,
        yend = y,
        color = matched
      ),
      show.legend = FALSE
    ) +
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
      size = .7 * cex_opt,
      linetype = 'longdash',
      show.legend = FALSE
    ) +
    scale_color_manual(values = c(`TRUE` = 'forestgreen', `FALSE` = 'red'))

  # Annotate the offset number
  plot_data = plot_data +
    ggrepel::geom_text_repel(
      data = expectation %>% filter(!matched),
      aes(
        x = x,
        y = y,
        label = round(offset, 2),
        color = matched
      ),
      nudge_x = 0,
      nudge_y = 0,
      size = 3 * cex_opt,
      show.legend = FALSE
    )

  # Add QC info
  QC = x$peaks_analysis$matches %>%
    dplyr::filter(karyotype == k) %>%
    dplyr::filter(row_number() == 1) %>%
    dplyr::pull(QC)

  # Function to a border to a plot
  qc_plot = function(x, QC)
  {
    if (is.na(QC))
      return(x)
    qc = ifelse(QC == "FAIL", "indianred3", 'forestgreen')

    x +
      theme(title = element_text(color = qc),
            panel.border = element_rect(colour = qc,
                                        fill = NA))
  }


  return(qc_plot(plot_data, QC))
}
