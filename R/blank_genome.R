blank_genome = function(chromosomes = paste0('chr', c(1:22, 'X', 'Y')))
{
  # GEt hg19 coordinates
  data('chr_coordinates_hg19', package = 'CNAqc')

  chr_coordinates_hg19 = chr_coordinates_hg19 %>% filter(chr %in% chromosomes)

  low = min(chr_coordinates_hg19$from)
  upp = max(chr_coordinates_hg19$to)

  ggplot(chr_coordinates_hg19) +
    my_ggplot_theme() +
    geom_rect(
      aes(
        xmin = centromerStart,
        xmax = centromerEnd,
        ymin = -Inf,
        ymax = Inf
      ),
      alpha = .3,
      colour = 'gainsboro'
    ) +
    geom_vline(xintercept = chr_coordinates_hg19$from,
               size = 0.3,
               colour = 'black') +
    geom_label(
      data = chr_coordinates_hg19,
      aes(
        x = chr_coordinates_hg19$from,
        y = -0.5,
        label = gsub('chr', '', chr_coordinates_hg19$chr)
      ),
      hjust = 0,
      colour = 'white',
      fill = 'black',
      size = 3
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
    labs(x = "Location",
         y = "Major/ minor allele") +
    ggpubr::rotate_y_text() +
    xlim(low, upp)
}

