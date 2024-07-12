blank_genome_multisample = function(x,
                                    ref = "GRCh38",
                                    chromosomes = paste0('chr', c(1:22, 'X', 'Y')),
                                    label_chr = -0.5,
                                    cex = 1) {
  reference_coordinates = get_reference(ref) %>%
    filter(chr %in% chromosomes)
  
  low = min(reference_coordinates$from)
  upp = max(reference_coordinates$to)
  sample_names = get_sample_name(x)
  
  reference_coordinates = map2_dfr(1:length(get_sample_name(x)), sample_names, ~ reference_coordinates %>% 
                                    mutate(sample_id = .y))
  
  p1 = ggplot2::ggplot(reference_coordinates, aes(x = from, y = sample_id)) +
    CNAqc:::my_ggplot_theme(cex = cex) +
    ggplot2::geom_linerange(
      ggplot2::aes(xmin = centromerStart,
                   xmax = centromerEnd,),
      size = .1,
      color = 'black',
      linetype = 8
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = centromerStart,
        xend = centromerEnd,
        y = 0,
        yend = Inf
      ),
      size = .1,
      color = 'black',
      linetype = 8
    )
  
  
  p1 = p1 + ggplot2::geom_rect(
    data = reference_coordinates,
    ggplot2::aes(
      xmin = from,
      xmax = from,
      ymin = 0,
      ymax = Inf
    ),
    alpha = 1,
    colour = 'grey',
  )
  
  p1 = p1 +
    # ggplot2::geom_hline(yintercept = 0,
    #                     size = 1,
    #                     colour = 'gainsboro') +
    # ggplot2::geom_hline(
    #   yintercept = 1,
    #   size = .3,
    #   colour = 'black',
    #   linetype = 'dashed'
    # ) +
    ggplot2::labs(x = "Chromosome",
                  y = "Samples") +
    # ggpubr::rotate_y_text() +
    # ggpubr::rotate_x_text() +
    # xlim(low, upp) +
    
    #set the chr names in the centromer positions.
    ggplot2::scale_x_continuous(
      breaks = c(0, reference_coordinates$centromerStart, upp),
      labels = c(
        "",
        gsub(
          pattern = 'chr',
          replacement = '',
          reference_coordinates$chr
        ),
        ""
      )
    )
  return(p1)
}
