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

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
plot_arm_fragmentation = function(x)
{
  # add a test: if there are any, etc


  # Results from fragmentation test
  fragments = x$arm_fragmentation$table %>% filter(significant)
  segments = CNAqc:::relative_to_absolute_coordinates(x$cna)

  # Plot the fragments table
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
      title = paste0("Overfragmentation (", nrow(fragments), ' arms)'),
      x = paste0('Segments >', 100 * x$arm_fragmentation$genome_percentage_cutoff,'% of arm'),
      y = paste0('Segments <', 100 * x$arm_fragmentation$genome_percentage_cutoff,'% of arm'),
      subtitle = paste0(
        "Bonferroni cut ", x$arm_fragmentation$table$Bonferroni_cutoff[1],
        ' (alpha = ', x$arm_fragmentation$alpha, ', ', x$arm_fragmentation$N_tests, ' tests).')
    ) +
    guides(size = guide_legend("-log(p)"), color = guide_legend("Fragmented")) +
    ggrepel::geom_text_repel(
      data = fragments,
      segment.size = 0.2,
      segment.color = "grey50",
      aes(label = paste0(gsub('chr', '', chr), arm), x = n_long, y = n_short),
      size = 3,
      color = 'forestgreen',
      xlim = 10
    )

  # Plot the jumps
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
      subtitle = paste0(
        "Maximum jump ", max(x$arm_fragmentation$table$jumps), '.'),
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

  # Zoom in each chromosome
  chr_to_plot = fragments$chr %>% unique

  plots_zoom = lapply(chr_to_plot,
         function(chr){
           # Add breaks for every segment
           segments_plot = CNAqc::plot_segments(x, chromosomes = chr)
           chr_breakpoints = segments %>% filter(chr == !!chr)

           segments_plot +
             labs(caption = NULL) +
             geom_vline(xintercept = chr_breakpoints$from, size = .2, linetype = 'dashed')
         })

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

  NP = ceiling(length(plots_zoom)/2)

  zoom_panel =
    ggarrange(
      plotlist = plots_zoom,
      nrow =  NP,
      ncol = 2
    )

  figure =
    ggarrange(
      top_panel,
      zoom_panel,
      nrow = 2,
      ncol = 1,
      heights = c(1.2, NP)
    )

  return(figure)
}
