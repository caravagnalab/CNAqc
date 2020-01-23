# x = detect_overfragmentation(x, genome_percentage_cutoff = .1)
# plot_arm_fragmentation(x)

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
      subtitle = paste0(
        "Bonferroni cut ", x$arm_fragmentation$table$Bonferroni_cutoff[1],
        ' (alpha = ', x$arm_fragmentation$alpha, ', ', x$arm_fragmentation$N_tests, ' tests).')
    ) +
    guides(size = guide_legend("-log(p)"), color = guide_legend("Fragmented")) +
    ggrepel::geom_text_repel(
      data = fragments,
      aes(label = paste0(gsub('chr', '', chr), arm), x = n_long, y = n_short),
      size = 3,
      color = 'forestgreen',
      xlim = 10
    )

  # Zoom in each chromosome
  chr_to_plot = fragments$chr %>% unique

  plots_zoom = lapply(chr_to_plot,
         function(chr){
           # Add breaks for every segment
           segments_plot = CNAqc::plot_segments(x, chromosomes = chr)
           chr_breakpoints = segments %>% filter(chr == !!chr)

           segments_plot +
             labs(caption = paste0("Chromsome zoom")) +
             geom_vline(xintercept = chr_breakpoints$from, size = .2, linetype = 'dashed')
         })

  NP = 1 + length(plots_zoom)
  M = ceiling(sqrt(NP))
  X = Y = M

  if(NP == 2) { X = 1; Y = 2 }

  figure =
    ggarrange(
      plotlist = append(list(fragments_table_plot), plots_zoom),
      nrow = X,
      ncol = Y
    )

  return(figure)
}
