blank_genome = function(chromosomes = paste0('chr', c(1:22, 'X', 'Y')), label_chr = -0.5)
{
  # GEt hg19 coordinates
  data('chr_coordinates_hg19', package = 'CNAqc')

  chr_coordinates_hg19 = chr_coordinates_hg19 %>% filter(chr %in% chromosomes)

  low = min(chr_coordinates_hg19$from)
  upp = max(chr_coordinates_hg19$to)

  pl = ggplot(chr_coordinates_hg19) +
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
      data = chr_coordinates_hg19,
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
      breaks = c(0, chr_coordinates_hg19$from, upp),
      labels = c("", gsub(pattern = 'chr', replacement = '', chr_coordinates_hg19$chr), "")
    )


  # if(!is.null(label_chr) & !is.na(label_chr))
  #   pl = pl +
  #   geom_label(
  #     data = chr_coordinates_hg19,
  #     aes(
  #       x = chr_coordinates_hg19$from,
  #       y = label_chr,
  #       label = gsub('chr', '', chr_coordinates_hg19$chr)
  #     ),
  #     hjust = 0,
  #     colour = 'white',
  #     fill = 'black',
  #     size = 3
  #   )


  return(pl)
}

