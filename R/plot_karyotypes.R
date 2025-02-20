#' Plot counts and numbers of clonal simple CNAs.
#'
#' @description Plot a barplot for the segments, reporting either their
#'  counts of the proportion of genome covered.
#'
#' @param x A CNAqc object.
#' @param type With \code{"percentage"}, the proportion of genome covered is returned.
#' With \code{"number"}, the segment counts. In all other cases an erorr is generated.
#' @param chromosomes The chromosome to use for this plot.
#'
#' @return A \code{ggplot2} plot.
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(mutations = example_dataset_CNAqc$mutations, cna = example_dataset_CNAqc$cna, purity = example_dataset_CNAqc$purity)
#'
#' # All chromosomes (default)
#' plot_karyotypes(x)
#'
#' # Specific chromosomes
#' plot_karyotypes(x, chromosomes = 'chr3')
#' plot_karyotypes(x, chromosomes = 'chr13')
plot_karyotypes = function(x,
                           type = "percentage",
                           chromosomes = paste0('chr', c(1:22, 'X', 'Y')))
{
  stopifnot(inherits(x, 'cnaqc'))

  # Get coordinates for used chromosomes
  genome_size = CNAqc:::get_reference(x$reference_genome, data = x$genomic_coordinates) %>%
    dplyr::filter(chr %in% chromosomes)

  segments = x$cna %>%
    dplyr::filter(chr %in% chromosomes) %>%
    dplyr::mutate(
      label = paste0(Major, ':', minor),
      call = ifelse(CCF == 1, "Clonal", "Subclonal"),
      size = to - from
    )

  colors = CNAqc:::get_karyotypes_colors('other')

  segments = segments %>%
    dplyr::mutate(label = ifelse(label %in% names(colors), label, 'other')) %>%
    select(size, label, call) %>%
    group_by(label, call) %>%
    summarise(size = sum(size))

  if(type == 'number')
  {
    pl = ggplot2::ggplot(segments,
                         ggplot2::aes(x = '', fill = label, y = size)) +
      my_ggplot_theme() +
      ggplot2::geom_bar(alpha = 1, color = 'white', size = .1, stat = 'identity') +
      ggplot2::scale_fill_manual(values = colors) +
      ggplot2::labs(y = "Count", x = 'Karyotypes', title = "Length of segments") +
      ggplot2::guides(fill = ggplot2::guide_legend('')) +
      ggplot2::facet_wrap(~call)

    return(pl)
  }
  else
  {
    # genome_size = sum(genome_size$length)

    segments = segments %>%
      mutate(p = size/sum(segments$size))

    # segments = segments %>%
    #   dplyr::filter(CCF == 1) %>%
    #   dplyr::mutate(
    #     percentage = (to - from)/genome_size
    #   ) %>%
    #   dplyr::group_by(label, call) %>%
    #   dplyr::summarise(percentage = sum(percentage))

    colors = get_karyotypes_colors(unique(segments$label))

    pl = ggplot2::ggplot(segments,
                         ggplot2::aes(x = '', y = p, fill = label)) +
      my_ggplot_theme() +
      ggplot2::ylim(0, 1) +
      ggplot2::geom_bar(stat = 'identity', alpha = 1, color = 'white', size = .1) +
      ggplot2::scale_fill_manual(values = colors) +
      ggplot2::labs(y = "Percentage", x = 'Karyotype', title = "Genome coverage") +
      ggplot2::guides(fill = ggplot2::guide_legend('')) +
      ggplot2::facet_wrap(~call)

    return(pl)
  }

  stop("Type = {percentage, number}.")
}



