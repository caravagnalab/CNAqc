#' Plot a barplot for the segments.
#'
#'  @description Plot a barplot for the segments, reporting either their
#'  counts of the proportion of genome covered.
#'
#' @param x An object of class \code{cnaqc}, created by the \code{init} function.
#' @param type With \code{"percentage"}, the proportion of genome covered is returned.
#' With \code{"number"}, the segment counts. In all other cases an erorr is generated.
#' @param chromosomes The chromosome to use for this plot.
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna,example_dataset_CNAqc$purity)
#'
#' plot_karyotypes(x)
#'
#' plot_karyotypes(x, chromosomes = 'chr3')
#' plot_karyotypes(x, chromosomes = 'chr13')
plot_karyotypes = function(x,
                           type = "percentage",
                           chromosomes = paste0('chr', c(1:22, 'X', 'Y')))
{
  stopifnot(inherits(x, 'cnaqc'))

  data('chr_coordinates_hg19', package = 'CNAqc')
  genome_size = chr_coordinates_hg19 %>%
    filter(chr %in% chromosomes) %>%
    rename(chr_size = length)

  segments = x$cna %>%
    filter(chr %in% chromosomes) %>%
    mutate(
      label = paste0(Major, ':', minor),
      call = ifelse(CCF == 1, "Clonal", "Subclonal")
    )

  if(type == 'number')
  {
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

    segments = segments %>%
      filter(CCF == 1) %>%
      mutate(
        percentage = (to - from)/genome_size
      ) %>%
      group_by(label, call) %>%
      summarise(percentage = sum(percentage))

    colors = get_karyotypes_colors(unique(segments$label))

    pl = ggplot(segments %>% mutate(K = ""),
           aes(x = K, y = percentage, fill = label)) +
      my_ggplot_theme() +
      ylim(0, 1) +
      geom_bar(stat = 'identity', alpha = 1, color = 'white', size = .1) +
      scale_fill_manual(values = colors) +
      labs(y = "Percentage", x = 'Karyotype', title = "Genome coverage") +
      guides(fill = guide_legend('')) +
      facet_wrap(~call)

    return(pl)
  }

  stop("Type = {percentage, number}.")
}



