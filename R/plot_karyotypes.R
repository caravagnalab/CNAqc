#' Title
#'
#' @param x
#' @param type
#' @param chromosomes
#'
#' @return
#' @export
#'
#' @examples
plot_karyotypes = function(x,
                           type = "percentage",
                           chromosomes = paste0('chr', c(1:22, 'X', 'Y')))
{
  stopifnot(inherits(x, 'cnaqc'))

  data('chr_coordinates_hg19', package = 'CNAqc')
  genome_size = chr_coordinates_hg19 %>%
    filter(chr %in% chromosomes) %>%
    rename(chr_size = length)

  if(type == 'number')
  {
    segments = x$cna %>%
      filter(chr %in% chromosomes) %>%
      mutate(
        label = paste0(Major, ':', minor),
        call = ifelse(CCF == 1, "Clonal", "Subclonal")
        )

    colors = get_karyotypes_colors(unique(segments$label))

    pl = ggplot(segments %>% mutate(K = ""),
                 aes(x = K, fill = label)) +
      my_ggplot_theme() +
      geom_bar(alpha = 1, color = 'white', size = .1) +
      scale_fill_manual(values = colors) +
      labs(y = "Count", x = 'Karyotypes', title = "Number of segments") +
      guides(fill = guide_legend('')) +
      facet_wrap(~call)

    return(pl)
  }
  else
  {
    genome_size = rev(chr_coordinates_hg19$to)[1]

    segments = x$cna %>%
      filter(
        CCF == 1,
        chr %in% chromosomes
      ) %>%
      mutate(
        label = paste0(Major, ':', minor),
        percentage = (to - from)/genome_size
      ) %>%
      group_by(label) %>%
      summarise(percentage = sum(percentage))

    colors = get_karyotypes_colors(unique(segments$label))

    pl = ggplot(segments %>% mutate(K = ""),
           aes(x = K, y = percentage, fill = label)) +
      my_ggplot_theme() +
      ylim(0, 1) +
      geom_bar(stat = 'identity', alpha = 1, color = 'white', size = .1) +
      scale_fill_manual(values = colors) +
      labs(y = "Percentage", x = 'Karyotype', title = "Genome coverage") +
      guides(fill = guide_legend(''))

    return(pl)
  }

  stop("Type = {percentage, number}.")
}



