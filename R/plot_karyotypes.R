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
#' x = init(example_dataset_CNAqc$mutations, example_dataset_CNAqc$cna,example_dataset_CNAqc$purity)
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

  # Get coordinates for used chromosomes
  genome_size = CNAqc:::get_reference(x$reference_genome) %>%
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
    pl = ggplot(segments,
                 aes(x = '', fill = label, y = size)) +
      my_ggplot_theme() +
      geom_bar(alpha = 1, color = 'white', size = .1, stat = 'identity') +
      scale_fill_manual(values = colors) +
      labs(y = "Count", x = 'Karyotypes', title = "Length of segments") +
      guides(fill = guide_legend('')) +
      facet_wrap(~call)

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

    pl = ggplot(segments,
           aes(x = '', y = p, fill = label)) +
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



