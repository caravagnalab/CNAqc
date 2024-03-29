#' Plot the arm level fragmentation test.
#'
#' @description
#'
#' This function produces a arm-level report for the fragmentation test, with:
#'
#' * a scatter of the counts per arm, with scaled the p-values;
#' * a jump statistics per arm, $J$,  the sum of the variation in total copy number
#'  profiles, evaluated among each pair of contiguous segments.
#'
#'  Significantly overfragmented arms with high $J$ have a “scattered” copy number
#'   profile. Those with low $J$ are more uniform, as they show little no copy
#'   number change, and can be possibly smoothed (see below).
#'
#' @param x A CNAqc object.
#' @param zoom Number of maximum zoom panels to show in the bottom of
#' the figure. By default 0.
#'
#' @return A \code{ggpubr} figure
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(mutations = example_dataset_CNAqc$mutations, cna = example_dataset_CNAqc$cna, purity = example_dataset_CNAqc$purity)
#'
#' x = detect_arm_overfragmentation(x, genome_percentage_cutoff = .2)
#' plot_arm_fragmentation(x)
plot_arm_fragmentation = function(x, zoom = 0)
{
  # add a test: if there are any, etc
  if (all(is.null(x$arm_fragmentation))) {
    cli::cli_alert_danger("Arm overfragmentation is not available for this object, compute it first.")
    return(CNAqc:::eplot())
  }

  # plot cex overriding
  cex_opt = getOption('CNAqc_cex', default = 1)

  # Results from fragmentation test
  fragments = x$arm_fragmentation$table %>% filter(significant)
  segments = CNAqc:::relative_to_absolute_coordinates(x, x$cna)

  ###### ###### ###### ###### ######
  # Plot the fragments table
  ###### ###### ###### ###### ######
  squaring = max(c(
    x$arm_fragmentation$table$n_long,
    x$arm_fragmentation$table$n_short
  ))
  if (squaring <= 10)
    squaring = 15

  # fragments_table_plot =
  #   ggplot(
  #   x$arm_fragmentation$table
  # ) +
  #   geom_point(aes(x = n_long, y = n_short, size = -log(p_value), color = significant)) +
  #   CNAqc:::my_ggplot_theme() +
  #   geom_abline(size = .3, color = 'red', linetype = 'dashed') +
  #   xlim(0, squaring) +
  #   ylim(0, squaring) +
  #   scale_color_manual(values = c(`FALSE` = 'indianred3', `TRUE` = 'forestgreen')) +
  #   labs(
  #     title = paste0(nrow(fragments), " overfragmented arms (Bonferroni alpha = ", x$arm_fragmentation$table$Bonferroni_cutoff[1], ')'),
  #     x = paste0('Long segments (>', 100 * x$arm_fragmentation$genome_percentage_cutoff,'%)'),
  #     y = paste0('Short segments (<', 100 * x$arm_fragmentation$genome_percentage_cutoff,'%)')
  #   ) +
  #   guides(size = guide_legend("-log(p)", nrow = 1), color = guide_legend("Fragmented")) +
  #   ggrepel::geom_text_repel(
  #     data = fragments,
  #     segment.size = 0.2,
  #     segment.color = "grey50",
  #     aes(label = paste0(gsub('chr', '', chr), arm), x = n_long, y = n_short),
  #     size = 3,
  #     color = 'forestgreen',
  #     xlim = 10
  #   )

  fragments_table_plot = ggplot2::ggplot(x$arm_fragmentation$table) +
    ggplot2::geom_point(ggplot2::aes(
      x = n_long,
      y = n_short,
      size = jumps,
      color = significant
    )) +
    ggplot2::scale_size_continuous(range = c(1, 3) * cex_opt)  +
    CNAqc:::my_ggplot_theme() +
    ggplot2::geom_abline(size = .3,
                         color = 'red',
                         linetype = 'dashed') +
    ggplot2::xlim(0, squaring) +
    ggplot2::ylim(0, squaring) +
    ggplot2::scale_color_manual(values = c(`FALSE` = 'orange', `TRUE` = 'steelblue')) +
    labs(
      title = bquote(
        .(nrow(fragments)) ~ "overfragmented arms at" ~ alpha ~ '=' ~ .(x$arm_fragmentation$alpha) ~
          'with FWER'
      ),
      x = paste0(
        'Long segments (>',
        100 * x$arm_fragmentation$genome_percentage_cutoff,
        '%)'
      ),
      y = paste0(
        'Short segments (<',
        100 * x$arm_fragmentation$genome_percentage_cutoff,
        '%)'
      ),
      caption = paste0(
        "Bonferroni adjusted alpha = ",
        x$arm_fragmentation$table$Bonferroni_cutoff[1],
        '~ ',
        sum(fragments$significant),
        ' significant tests.'
      )
    ) +
    ggplot2::guides(size = ggplot2::guide_legend("Jump", nrow = 1), color = 'none') +
    ggrepel::geom_text_repel(
      data = fragments,
      segment.size = 0.2,
      segment.color = "grey50",
      aes(
        label = paste0(gsub('chr', '', chr), arm),
        x = n_long,
        y = n_short
      ),
      size = 4.5 * cex_opt,
      color = 'steelblue',
      xlim = 10
    ) +
    ggplot2::theme(plot.caption = ggplot2::element_text(color = ifelse(
      sum(fragments$significant) > 0, "steelblue",  "orange"
    )))

  # ###### ###### ###### ###### ######
  # # Plot the jumps
  # ###### ###### ###### ###### ######
  # MJ = max(x$arm_fragmentation$table$jumps) - min(x$arm_fragmentation$table$jumps)
  #
  # J_order_x = lapply(
  #   paste0('chr', c(1:22, 'X', 'Y')),
  #   function(x) paste0(x, c('p', 'q'))
  # )
  # J_order_x = Reduce(c, J_order_x)
  #
  # plots_jumps = ggplot(
  #   x$arm_fragmentation$table
  # ) +
  #   geom_bar(aes(x = paste0(chr, arm), y = jumps, fill = significant), stat = 'identity') +
  #   CNAqc:::my_ggplot_theme() +
  #   scale_fill_manual(values = c(`FALSE` = 'indianred3', `TRUE` = 'forestgreen')) +
  #   labs(
  #     title = paste0("Segments jumps"),
  #     x = ''
  #   ) +
  #   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #   scale_x_discrete(limits = J_order_x) +
  #   coord_cartesian(clip = 'off')
  #   # guides(size = guide_legend("-log(p)"), color = guide_legend("Fragmented")) +
  #
  # if(MJ > 100)
  #   plots_jumps = plots_jumps +
  #     scale_y_log10() +
  #     labs(y = 'jumps (logscale)')
  #
  # ###### ###### ###### ###### ######
  # # Zoom plot if required
  # ###### ###### ###### ###### ######
  # plots_zoom = NULL

  if (zoom > 0)
  {
    # Zoom in each chromosome
    chr_to_plot = fragments$chr %>% unique

    if (length(chr_to_plot) > zoom) {
      cli::cli_alert_danger("Increare zoom to see more chromosomes.")
      chr_to_plot = chr_to_plot[1:zoom]
    }


    max_Y_height = 6

    plots_zoom = lapply(chr_to_plot,
                        function(chr) {
                          # Add breaks for every segment
                          segments_plot = CNAqc::plot_segments(
                            x,
                            chromosomes = chr,
                            circular = T,
                            max_Y_height = max_Y_height
                          )
                          chr_breakpoints = segments %>%
                            filter(chr == !!chr) %>%
                            mutate(y = Inf,
                                   outern = Major > max_Y_height)

                          j = x$arm_fragmentation$table %>% filter(chr == !!chr) %>% pull(jumps)

                          segments_plot +
                            ggplot2::labs(
                              caption = NULL,
                              title = chr,
                              subtitle = paste0(
                                'B = ',
                                nrow(chr_breakpoints),
                                ', J = ',
                                paste(j, collapse = ' + ')
                              )
                            ) +
                            ggplot2::geom_vline(
                              xintercept = chr_breakpoints %>%
                                filter(outern) %>%
                                pull(from),
                              size = .2,
                              linetype = 'dashed',
                              color = 'orange'
                            ) +
                            ggplot2::theme(plot.margin = ggplot2::margin(0.2, 0, 0, 0, "cm")) +
                            ggplot2::theme_void()
                        })

  }


  ###### ###### ###### ###### ######
  # Figure assembly
  ###### ###### ###### ###### ######

  top_panel = fragments_table_plot
  # top_panel = cowplot::plot_grid(
  #   fragments_table_plot + theme(legend.position="none"),
  #   plots_jumps + theme(legend.position="none"),
  #   nrow = 1,
  #   ncol = 2,
  #   align = 'h',
  #   rel_widths = c(.7, 1)
  # )

  # extract a legend that is laid out horizontally
  # legend_b <- get_legend(
  #   fragments_table_plot +
  #     theme(legend.position = "bottom")
  # )

  # add the legend underneath the row we made earlier. Give it 10%
  # of the height of one plot (via rel_heights).
  # top_panel = cowplot::plot_grid(top_panel, legend_b, ncol = 1, rel_heights = c(1, .1))
  figure = top_panel

  if (zoom > 0)
  {
    NP = ceiling(length(plots_zoom) / 4)

    zoom_panel =
      ggpubr::ggarrange(plotlist = plots_zoom,
                        nrow =  NP,
                        ncol = 4)

    figure =
      ggpubr::ggarrange(
        top_panel,
        zoom_panel,
        nrow = 2,
        ncol = 1,
        heights = c(1.2, NP)
      )
  }

  return(figure)
}

# plot_arm_fragmentation3d = function(x)
# {
#   # add a test: if there are any, etc
#   if(all(is.null(x$arm_fragmentation))) {
#     cli::cli_alert_danger("Arm overfragmentation is not available for this object, compute it first.")
#     return(CNAqc:::eplot())
#   }
#
#   require(plotly)
#
#   p <- plot_ly(
#     x$arm_fragmentation$table,
#     x = ~n_long, y = ~n_short, z = ~jumps, color = ~significant, colors = c(`TRUE` = '#BF382A', `FALSE`='#0C4B8E')) %>%
#     add_markers() %>%
#     layout(scene = list(xaxis = list(title = 'Long'),
#                         yaxis = list(title = 'Short'),
#                         zaxis = list(title = 'Jumps')))
#
#   p
# }
