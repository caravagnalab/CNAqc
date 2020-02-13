blank_genome = function(ref = "GRCh38", chromosomes = paste0('chr', c(1:22, 'X', 'Y')), label_chr = -0.5)
{
  reference_coordinates = get_reference(ref) %>% 
    filter(chr %in% chromosomes)

  low = min(reference_coordinates$from)
  upp = max(reference_coordinates$to)

  pl = ggplot(reference_coordinates) +
    CNAqc:::my_ggplot_theme() +
    geom_rect(
      aes(
        xmin = centromerStart,
        xmax = centromerEnd,
        ymin = 0,
        ymax = Inf
      ),
      alpha = .3,
      colour = 'gainsboro'
    ) +
    geom_segment(
      data = reference_coordinates,
      aes(
        x = from,
        xend = from,
        y = 0,
        yend = Inf
      ),
      size = .1,
      color = 'black',
      linetype = 8
    ) +
    geom_hline(yintercept = 0,
               size = 1,
               colour = 'gainsboro') +
    geom_hline(
      yintercept = 1,
      size = .3,
      colour = 'black',
      linetype = 'dashed'
    ) +
    labs(x = "Chromosome",
         y = "Major/ minor allele") +
    ggpubr::rotate_y_text() +
    # ggpubr::rotate_x_text() +
    # xlim(low, upp) +
    scale_x_continuous(
      breaks = c(0, reference_coordinates$from, upp),
      labels = c("", gsub(pattern = 'chr', replacement = '', reference_coordinates$chr), "")
    )


  # if(!is.null(label_chr) & !is.na(label_chr))
  #   pl = pl +
  #   geom_label(
  #     data = reference_coordinates,
  #     aes(
  #       x = reference_coordinates$from,
  #       y = label_chr,
  #       label = gsub('chr', '', reference_coordinates$chr)
  #     ),
  #     hjust = 0,
  #     colour = 'white',
  #     fill = 'black',
  #     size = 3
  #   )


  return(pl)
}

