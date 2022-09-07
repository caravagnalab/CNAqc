#' Plot the results of peak analysis.
#'
#'  @description Results from \code{analyze_peaks} can be visualised with this
#'  function, which arranges the plots of each karyotype in a figure via \code{ggpubr}.
#'  Each karyotype shows the data, the estimated density, the peaks (selected and
#'  discarded), and the fit with shaded matching area.
#'
#' @param x An object of class \code{cnaqc}, where function \code{analyze_peaks} has
#' been computed.
#' @param empty_plot If data for one karyotype is missing, an empty plot is returned.
#' Otherwise the plot is not returned (NULL is forwarded).
#' @param assembly_plot If \code{TRUE}, a unique figure is returned with all the
#' plots assembled. Otherwise a list of plot is returned.
#' @param what What karyotypes should be plot. Value `common` or `simple` refers to clonal karyotypes used
#' for sample-level QC. `general` or `complex` for all the others. `subclonal` is for subclonal segments.
#'
#' @return A \code{ggpubr} object for an assembled figure.
#' @export
#'
#' @import ggpubr
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(example_dataset_CNAqc$mutations, example_dataset_CNAqc$cna,example_dataset_CNAqc$purity)
#'
#' x = analyze_peaks(x)
#' plot_peaks_analysis(x)
plot_peaks_analysis = function(x,
                               empty_plot = TRUE,
                               assembly_plot = TRUE,
                               what = "simple")
{
  stopifnot(inherits(x, "cnaqc"))

  if (what %in% c('simple', 'common'))
  {
    with_peaks = all(!is.null(x$peaks_analysis))
    if (!with_peaks) {
      warning("Input does not have peaks, see ?peaks_analysis to run peaks analysis.")
      return(CNAqc:::eplot())
    }

    karyotypes = x$peaks_analysis$fits %>% names

    order_karyotypes = c('1:0', '1:1', '2:0', '2:1', '2:2')
    # karyotypes = order_karyotypes[order_karyotypes %in% karyotypes]
    karyotypes = order_karyotypes

    # Plot each one of the fits
    plots = lapply(karyotypes, function(k) {
      if (all(is.null(x$peaks_analysis$fits[[k]]$matching)))
      {
        if (empty_plot)
          return(CNAqc:::eplot())
        else
          return(NULL)
      }

      return(suppressWarnings(suppressMessages(plot_peaks_fit(x, k))))
    })

    plots = plots[!sapply(plots, is.null)]
    if (length(plots) == 0) {
      cli::cli_alert_warning("Nothing to plot")
      return(CNAqc:::eplot())
    }

    # Overall QC
    qc = ifelse(x$peaks_analysis$QC == 'PASS', 'forestgreen', 'indianred3')

    # Plots assembly
    if (assembly_plot)
      plots = suppressWarnings(suppressMessages(
        ggpubr::ggarrange(
          plotlist = plots,
          nrow = 1,
          ncol = length(plots)
        ) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(color = qc),
            panel.border = ggplot2::element_rect(colour = qc,
                                                 fill = NA)
          )
      ))

    return(plots)
  }

  if (what %in% c('complex', 'general'))
  {
    with_peaks = all(!is.null(x$peaks_analysis$general))

    if (!with_peaks) {
      warning("Input does not have peaks, see ?peaks_analysis to run peaks analysis.")
      return(eplot())
    }

    return(plot_peaks_fit_general(x))
  }

  if (what == 'subclonal')
  {
    with_peaks = all(!is.null(x$peaks_analysis$subclonal))

    if (!with_peaks) {
      warning("Input does not have peaks, see ?peaks_analysis to run peaks analysis.")
      return(CNAqc:::eplot())
    }

    pl = plot_peaks_fit_subclonal(x)
    if (assembly_plot)
      pl = ggpubr::ggarrange(
        plotlist = pl,
        ncol = 1,
        nrow = length(pl),
        common.legend = TRUE,
        legend = 'bottom'
      )

    return(pl)
  }

}

# Plot a single run results with the standard karyotypes model
plot_peaks_fit = function(x, k)
{
  matching = x$peaks_analysis$matching_strategy

  ranges = x$peaks_analysis$matches %>%
    dplyr::filter(karyotype == k) %>%
    pull(epsilon)

  # Required input values
  mutations = x$mutations %>%
    dplyr::filter(karyotype == k) %>%
    dplyr::mutate(karyotype = paste0(karyotype, " (", matching, ")"))

  den = x$peaks_analysis$fits[[k]]$density
  expectation = x$peaks_analysis$fits[[k]]$matching %>%
    dplyr::mutate(karyotype = paste0(karyotype, " (", matching, ")"))

  xy_peaks = x$peaks_analysis$fits[[k]]$xy_peaks
  purity_error = x$peaks_analysis$purity_error

  # linear combination of the weight, split by number of peaks to match
  weight = x$peaks_analysis$matches %>%
    dplyr::filter(karyotype == k) %>%
    dplyr::pull(weight) %>%
    sum

  # Plots cex for anything that is not the main theme
  cex_opt = getOption('CNAqc_cex', default = 1)

  # Add QC info
  QC = x$peaks_analysis$matches %>%
    dplyr::filter(karyotype == k) %>%
    dplyr::filter(row_number() == 1) %>%
    dplyr::pull(QC)

  qc_color = ifelse(QC == "FAIL", "indianred3", 'forestgreen')


  title = bquote(bold(.(k)) ~
                   .(paste0(
                     ' (n = ', nrow(mutations), ', ', round(weight * 100, 1),  '%)'
                   )))

  # Plot the data
  plot_data =
    ggplot2::ggplot(data = mutations, aes(VAF)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ..density..),
                            binwidth = 0.01,
                            alpha = .3) +
    ggplot2::geom_line(
      data = data.frame(x = den$x, y = den$y),
      ggplot2::aes(x = x, y = y),
      size = .3,
      color = 'black'
    ) +
    CNAqc:::my_ggplot_theme() +
    ggplot2::labs(title = title,
                  y = 'KDE',
                  x = "VAF") +
    ggplot2::theme(legend.position = 'bottom')  +
    ggplot2::xlim(-0.01, 1.01) +
    ggplot2::facet_wrap( ~ karyotype) +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = qc_color))

  # Add points for peaks to plot
  plot_data = plot_data +
    ggplot2::geom_point(
      data = xy_peaks,
      ggplot2::aes(
        x = x,
        y = y,
        shape = discarded,
        size = counts_per_bin
      ),
      show.legend = FALSE
    ) +
    ggplot2::scale_shape_manual(values = c(`TRUE` = 1, `FALSE` = 16)) +
    ggplot2::scale_size(range = c(1, 3) * cex_opt)

  # Add expectation peaks, and matching colors
  plot_data = plot_data +
    ggplot2::geom_point(
      data = expectation,
      ggplot2::aes(x = x, y = y, color = matched),
      size = 2 * cex_opt,
      shape = 4,
      show.legend = FALSE
    ) +
    ggplot2::annotate(
      geom = 'rect',
      xmin = expectation$x - expectation$VAF_tolerance,
      xmax = expectation$x + expectation$VAF_tolerance,
      ymin = 0,
      ymax = Inf,
      color = NA,
      alpha = .4,
      fill = 'purple4'
    ) +
    ggplot2::geom_segment(
      data = expectation,
      ggplot2::aes(
        x = x,
        y = y,
        xend = peak,
        yend = y,
        color = matched
      ),
      show.legend = FALSE
    ) +
    ggplot2::annotate(
      geom = 'rect',
      xmin = expectation$peak - expectation$epsilon,
      xmax = expectation$peak + expectation$epsilon,
      ymin = 0,
      ymax = Inf,
      color = NA,
      alpha = .4,
      fill = 'steelblue'
    ) +
    ggplot2::geom_vline(
      data = expectation,
      ggplot2::aes(xintercept = peak, color = matched),
      size = .7 * cex_opt,
      linetype = 'longdash',
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(values = c(`TRUE` = 'forestgreen', `FALSE` = 'red'))

  # Annotate the offset number
  plot_data = plot_data +
    ggrepel::geom_text_repel(
      data = expectation %>% filter(!matched),
      ggplot2::aes(
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

  # # Function to a border to a plot
  # qc_plot = function(x, QC)
  # {
  #   if (is.na(QC))
  #     return(x)
  #   qc = ifelse(QC == "FAIL", "indianred3", 'forestgreen')
  #
  #   x +
  #     theme(plot.title = element_text(color = qc))
  # }
  #
  #
  # return(qc_plot(plot_data, QC))

  return(plot_data)

}

# Plot general peaks analysis
plot_peaks_fit_general = function(x)
{
  add_counts = function(w) {
    w$karyotype = paste0(w$karyotype, ' (n = ', x$n_karyotype[w$karyotype], ')')
    w %>% as_tibble()
  }

  analysis = x$peaks_analysis$general$analysis
  n_min = x$peaks_analysis$general$params$n_min
  epsilon =  x$peaks_analysis$general$params$epsilon
  expected_peaks = x$peaks_analysis$general$expected_peaks %>% add_counts()
  n_bootstrap = x$peaks_analysis$general$params$n_bootstrap
  data_peaks = x$peaks_analysis$general$data_peaks %>% add_counts()
  data_densities = x$peaks_analysis$general$data_densities  %>% add_counts()

  # plotting
  x$mutations %>%
    filter(karyotype %in% analysis) %>%
    add_counts() %>%
    ggplot2::ggplot(aes(VAF)) +
    ggplot2::geom_histogram(aes(y = ..density..), binwidth = 0.01, fill = 'gray') +
    ggplot2::facet_wrap( ~ karyotype, scales = 'free_y') +
    CNAqc:::my_ggplot_theme() +
    ggplot2::geom_line(data = data_densities,
                       ggplot2::aes(x = x, y = y),
                       inherit.aes = FALSE) +
    ggplot2::geom_vline(
      data = expected_peaks,
      ggplot2::aes(xintercept = peak, color = matched),
      linetype = 'dashed',
      show.legend = FALSE
    ) +
    ggplot2::geom_point(data = data_peaks,
                        ggplot2::aes(x = x, y = y),
                        inherit.aes = FALSE) +
    ggplot2::scale_color_manual(values = c(`FALSE` = 'indianred3', `TRUE` = 'forestgreen')) +
    ggplot2::scale_fill_manual(values = c(`FALSE` = 'indianred3', `TRUE` = 'forestgreen')) +
    ggplot2::geom_rect(
      data = data.frame(
        xmin = expected_peaks$peak - epsilon,
        xmax = expected_peaks$peak + epsilon,
        ymin = 0,
        ymax = Inf,
        matched = expected_peaks$matched,
        karyotype = expected_peaks$karyotype
      ),
      ggplot2::aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax,
        fill = matched
      ),
      inherit.aes = FALSE,
      alpha = .3
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend("Matched peak", override.aes = ggplot2::aes(alpha = 1))) +
    ggplot2::labs(
      title = bquote(
        "Generalised peak detection ("
        * n['min'] * ' > ' * .(n_min) * ', ' * epsilon * ' = ' * .(epsilon *
                                                                     100) * '%)'
      ),
      caption = bquote(n['nbootstrap'] * ' = ' * .(n_bootstrap))
    )
}

# Plot subclonal peaks analysis
plot_peaks_fit_subclonal = function(x)
{
  expected_peaks = x$peaks_analysis$subclonal$expected_peaks
  data_peaks = x$peaks_analysis$subclonal$data_peaks
  data_densities = x$peaks_analysis$subclonal$data_densities
  decision_table = x$peaks_analysis$subclonal$summary
  n_min = x$peaks_analysis$subclonal$params$n_min
  n_bootstrap = x$peaks_analysis$subclonal$params$n_bootstrap
  subclonal_mutations = x$peaks_analysis$subclonal$mutations
  epsilon =  x$peaks_analysis$subclonal$params$epsilon

  plot_model_id = function(segment_id)
  {
    this_model_peaks = expected_peaks %>% filter(segment_id == !!segment_id)
    this_model_ids = this_model_peaks$model_id %>% unique

    rank_models = decision_table %>%
      filter(segment_id == !!segment_id) %>%
      arrange(desc(prop)) %>%
      pull(model_id)

    rank_models = c(rank_models, setdiff(this_model_ids, rank_models))

    which_best = decision_table %>%
      filter(segment_id == !!segment_id) %>%
      arrange(desc(prop))
    which_best = which_best %>% filter(prop == which_best$prop[1]) %>% pull(model_id)

    strip_colors = rep("gray", rank_models %>% length())
    names(strip_colors) = rank_models
    strip_colors[which_best] = "goldenrod3"

    rep_muts = lapply(this_model_ids, function(x) {
      subclonal_mutations %>%
        filter(segment_id == !!segment_id) %>%
        mutate(model_id = x,
               model = ifelse(grepl('->', x), "linear", 'branching'))
    }) %>% Reduce(f = bind_rows)

    # rep_muts$model_id

    this_title = decision_table %>%
      filter(segment_id == !!segment_id) %>%
      filter(row_number() == 1) %>%
      select(segment_id, size, clones) %>% unlist() %>% paste(collapse = ' ')

    rep_muts %>%
      ggplot2::ggplot()  +
      ggplot2::geom_histogram(aes(x = VAF, y = ..density..),
                              binwidth = 0.01,
                              fill = 'gray') +
      ggplot2::xlim(-0.1, 1.1) +
      CNAqc:::my_ggplot_theme() +
      ggplot2::geom_rect(
        data = this_model_peaks,
        ggplot2::aes(
          xmin = peak - epsilon,
          xmax = peak + epsilon,
          fill = matched
        ),
        inherit.aes = FALSE,
        alpha = .3,
        ymin = 0,
        ymax = Inf
      ) +
      ggplot2::scale_color_manual(values = c(`FALSE` = 'indianred3', `TRUE` = 'forestgreen')) +
      ggplot2::scale_fill_manual(values = c(`FALSE` = 'indianred3', `TRUE` = 'forestgreen')) +
      ggplot2::geom_line(
        data = data_densities %>% filter(segment_id == !!segment_id),
        ggplot2::aes(x = x, y = y),
        inherit.aes = FALSE
      ) +
      ggplot2::geom_point(data = data_peaks %>% filter(segment_id == !!segment_id),
                          ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_vline(
        data = expected_peaks %>% filter(segment_id == !!segment_id),
        ggplot2::aes(
          xintercept = peak,
          linetype = role,
          color = matched
        )
      ) +
      ggplot2::facet_wrap( ~ factor(model_id, levels = rank_models)) +
      ggplot2::labs(title = this_title)


    # if(subclonal_mutations$segment_id %>% unique() %>% length() > 5)
    #   pl = pl + facet_wrap(segment_id~model, scales = 'free_y')
    # else
    #   pl = pl + facet_grid(segment_id~model, scales = 'free_y') +
    #   theme(strip.text.y.right = element_text(angle = 0))
  }

  decision_table$segment_id %>% unique() %>% lapply(plot_model_id)

  # what_model= expected_peaks$model_id %>% unique()

  # plotting
  # s1 = subclonal_mutations %>% mutate(model = 'branching')
  # s2 = subclonal_mutations %>% mutate(model = 'linear')
  #
  # pl = bind_rows(s1,s2) %>%
  #   ggplot()  +
  #   geom_histogram(aes(x = VAF, y = ..density..), binwidth = 0.01, fill = 'gray') +
  #   xlim(-0.1, 1.1) +
  #   CNAqc:::my_ggplot_theme()
  #
  # if(subclonal_mutations$segment_id %>% unique() %>% length() > 5)
  #   pl = pl + facet_wrap(segment_id~model, scales = 'free_y')
  # else
  #   pl = pl + facet_grid(segment_id~model, scales = 'free_y') +
  #   theme(strip.text.y.right = element_text(angle = 0))
  #
  # pl +
  #   geom_rect(
  #     data = data.frame(
  #       xmin = expected_peaks$peak - epsilon,
  #       xmax = expected_peaks$peak + epsilon,
  #       ymin = 0,
  #       ymax = Inf,
  #       model = expected_peaks$model,
  #       matched = expected_peaks$matched,
  #       segment_id = expected_peaks$segment_id
  #     ),
  #     aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = matched),
  #     inherit.aes = FALSE,
  #     alpha = .3
  #   ) +
  #   scale_color_manual(
  #     values = c(`FALSE` = 'indianred3', `TRUE` = 'forestgreen')
  #   ) +
  #   scale_fill_manual(
  #     values = c(`FALSE` = 'indianred3', `TRUE` = 'forestgreen')
  #   ) +
  #   geom_line(data = data_densities,
  #             aes(x = x, y = y),
  #             inherit.aes = FALSE) +
  #   geom_point(data = data_peaks,
  #              aes(x = x, y = y)) +
  #   geom_vline(data = expected_peaks,
  #              aes(xintercept = peak, linetype = role, color = matched))
}
