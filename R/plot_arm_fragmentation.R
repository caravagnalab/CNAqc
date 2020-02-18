# data('example_dataset_CNAqc')
# x = CNAqc::init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna, example_dataset_CNAqc$purity)
# x = detect_overfragmentation(x, genome_percentage_cutoff = .1)
# ggsave(plot_arm_fragmentation(x),filename = 'a.pdf', width = 10, height = 10)

#
# x = evoverse.datasets::MSEQ_CRC_ADENOCARCINOMA_SET7
# muts = x$mutations %>% select(chr, from, to, ref, alt, starts_with('Set7_55'))
# colnames(muts)[6:8] = c('VAF', 'DP', 'NV')
# cna = x$CNA %>% select(chr, from, to, starts_with('Set7_55'))
# colnames(cna)[4:5] = c('minor', 'Major')
# x = CNAqc::init(muts, cna, .88)
# CNAqc::plot_segments(x)
# x = detect_overfragmentation(x, genome_percentage_cutoff = .2)
# plot_arm_fragmentation(x)

#' Plot patterns of arm level fragmentation.
#'
#' @description
#'
#' The function produces a multi-panel figure. The top left
#' panel shows a scatter of the segments' length per arm, and
#' their jump value. The bottom panel uses a circular layout
#' to show arms, and can be activated setting \code{zoom > 0} 
#' (default is 0). In this case dashed lines outgoing the centre
#' reprent breakpoints not shown in the plot (i.e., with Major > 5).
#'
#' @param x A `CNAqc` object.
#' @param zoom Number of maximum zoom panels to show in the bottom of
#' the figure.
#'
#' @return
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna,example_dataset_CNAqc$purity)
#'
#' x = detect_arm_overfragmentation(x, genome_percentage_cutoff = .2)
#' plot_arm_fragmentation(x)
plot_arm_fragmentation = function(x, zoom = 0)
{
  # add a test: if there are any, etc
  if(all(is.null(x$arm_fragmentation))) {
    cli::cli_alert_danger("Arm overfragmentation is not available for this object, compute it first.")
    return(ggplot() + geom_blank())
  }

  # Results from fragmentation test
  fragments = x$arm_fragmentation$table %>% filter(significant)
  segments = CNAqc:::relative_to_absolute_coordinates(x, x$cna)

  ###### ###### ###### ###### ###### 
  # Plot the fragments table
  ###### ###### ###### ###### ###### 
  squaring = max(c(x$arm_fragmentation$table$n_long, x$arm_fragmentation$table$n_short))
  if(squaring <= 10) squaring = 15

  fragments_table_plot =
    ggplot(
    x$arm_fragmentation$table
  ) +
    geom_point(aes(x = n_long, y = n_short, size = -log(p_value), color = significant)) +
    CNAqc:::my_ggplot_theme() +
    geom_abline(size = .3, color = 'red', linetype = 'dashed') +
    xlim(0, squaring) +
    ylim(0, squaring) +
    scale_color_manual(values = c(`FALSE` = 'indianred3', `TRUE` = 'forestgreen')) +
    labs(
      title = paste0(nrow(fragments), " overfragmented arms (Bonferroni alpha = ", x$arm_fragmentation$table$Bonferroni_cutoff[1], ')'),
      x = paste0('Long segments (>', 100 * x$arm_fragmentation$genome_percentage_cutoff,'%)'),
      y = paste0('Short segments (<', 100 * x$arm_fragmentation$genome_percentage_cutoff,'%)')
    ) +
    guides(size = guide_legend("-log(p)", nrow = 1), color = guide_legend("Fragmented")) +
    ggrepel::geom_text_repel(
      data = fragments,
      segment.size = 0.2,
      segment.color = "grey50",
      aes(label = paste0(gsub('chr', '', chr), arm), x = n_long, y = n_short),
      size = 3,
      color = 'forestgreen',
      xlim = 10
    )

  ###### ###### ###### ###### ###### 
  # Plot the jumps
  ###### ###### ###### ###### ###### 
  MJ = max(x$arm_fragmentation$table$jumps) - min(x$arm_fragmentation$table$jumps)

  J_order_x = lapply(
    paste0('chr', c(1:22, 'X', 'Y')),
    function(x) paste0(x, c('p', 'q'))
  )
  J_order_x = Reduce(c, J_order_x)

  plots_jumps = ggplot(
    x$arm_fragmentation$table
  ) +
    geom_bar(aes(x = paste0(chr, arm), y = jumps, fill = significant), stat = 'identity') +
    CNAqc:::my_ggplot_theme() +
    scale_fill_manual(values = c(`FALSE` = 'indianred3', `TRUE` = 'forestgreen')) +
    labs(
      title = paste0("Segments jumps"),
      x = ''
    ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_x_discrete(limits = J_order_x) +
    coord_cartesian(clip = 'off')
    # guides(size = guide_legend("-log(p)"), color = guide_legend("Fragmented")) +

  if(MJ > 100)
    plots_jumps = plots_jumps +
      scale_y_log10() +
      labs(y = 'jumps (logscale)')

  ###### ###### ###### ###### ###### 
  # Zoom plot if required
  ###### ###### ###### ###### ###### 
  plots_zoom = NULL
  
  if(zoom > 0)
  {
    # Zoom in each chromosome
    chr_to_plot = fragments$chr %>% unique
    
    if(length(chr_to_plot) > 9) {
      cli::cli_alert_danger(
        "Arm overfragmentation involves >9 arms, but only 9 (3x3) will be plot."
      )
      chr_to_plot = chr_to_plot[1:9]
    }
    
    
    max_Y_height = 6
    
    plots_zoom = lapply(chr_to_plot,
                        function(chr){
                          # Add breaks for every segment
                          segments_plot = CNAqc::plot_segments(x, chromosomes = chr, circular = T, max_Y_height = max_Y_height)
                          chr_breakpoints = segments %>%
                            filter(chr == !!chr) %>%
                            mutate(
                              y = Inf,
                              outern = Major > max_Y_height
                            )
                          
                          j = x$arm_fragmentation$table %>% filter(chr == !!chr) %>% pull(jumps)
                          
                          segments_plot +
                            labs(
                              caption = NULL,
                              title = chr,
                              subtitle = paste0('B = ', nrow(chr_breakpoints), ', J = ', paste(j, collapse = ' + '))
                            ) +
                            geom_vline(xintercept = chr_breakpoints %>%
                                         filter(outern) %>%
                                         pull(from), size = .2, linetype = 'dashed', color = 'orange') +
                            theme(plot.margin = margin(0.2, 0, 0, 0, "cm")) +
                            theme_void()
                        })
    
  }


  ###### ###### ###### ###### ###### 
  # Figure assembly
  ###### ###### ###### ###### ###### 
  
  top_panel = cowplot::plot_grid(
    fragments_table_plot + theme(legend.position="none"),
    plots_jumps + theme(legend.position="none"),
    nrow = 1,
    ncol = 2,
    align = 'h',
    rel_widths = c(.7, 1)
  )

  # extract a legend that is laid out horizontally
  legend_b <- get_legend(
    fragments_table_plot +
      theme(legend.position = "bottom")
  )

  # add the legend underneath the row we made earlier. Give it 10%
  # of the height of one plot (via rel_heights).
  top_panel = cowplot::plot_grid(top_panel, legend_b, ncol = 1, rel_heights = c(1, .1))
  figure = top_panel
  
  if(zoom > 0)
  {
    NP = ceiling(length(plots_zoom)/4)
  
    zoom_panel =
      ggarrange(
        plotlist = plots_zoom,
        nrow =  NP,
        ncol = 4
      )
  
    figure =
      ggarrange(
        top_panel,
        zoom_panel,
        nrow = 2,
        ncol = 1,
        heights = c(1.2, NP)
      )
  }
  
  return(figure)
}

plot_arm_fragmentation3d = function(x)
{
  # add a test: if there are any, etc
  if(all(is.null(x$arm_fragmentation))) {
    cli::cli_alert_danger("Arm overfragmentation is not available for this object, compute it first.")
    return(ggplot() + geom_blank())
  }

  require(plotly)

  p <- plot_ly(
    x$arm_fragmentation$table,
    x = ~n_long, y = ~n_short, z = ~jumps, color = ~significant, colors = c(`TRUE` = '#BF382A', `FALSE`='#0C4B8E')) %>%
    add_markers() %>%
    layout(scene = list(xaxis = list(title = 'Long'),
                        yaxis = list(title = 'Short'),
                        zaxis = list(title = 'Jumps')))

  p
}



