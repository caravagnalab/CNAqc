#' A circular plot for simple clonal CNAs.
#'
#' @description
#'
#' An icon is a small circular plot with no legend or axis, which
#' shows simple clonal CNAs.
#'
#' @param x A CNAqc object.
#'
#' @return A `ggplot2` plot.
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(mutations = example_dataset_CNAqc$mutations, cna = example_dataset_CNAqc$cna, purity = example_dataset_CNAqc$purity)
#'
#' plot_icon_CNA(x)
plot_icon_CNA = function(x)
{
  stopifnot(inherits(x, 'cnaqc'))

  L = x

  KARYO_colors =
    c(
      `0:0` = 'darkblue',
      `1:0` = 'steelblue',
      `1:1` = ggplot2::alpha('forestgreen', .8),
      `2:0` = 'turquoise4',
      `2:1` = ggplot2::alpha('orange', .8),
      `2:2` = 'indianred3',
      `other` = 'darkred'
    )


  # Extract calls, and flatten them for plotting
  calls = L$cna %>%
    dplyr::mutate(
      label = paste(Major, minor, sep = ':'),
      CN = minor + Major
    ) %>%
    dplyr::select(chr, from, to, label, CN)


  calls = CNAqc:::relative_to_absolute_coordinates(x, calls)

  calls_flat = calls %>%
    dplyr::mutate(
          label = ifelse(label %in% names(KARYO_colors), label, 'other')
        )

  chromosomes = calls_flat$chr %>% unique

  # Get coordinates for used chromosomes
  reference_genome = CNAqc:::get_reference(x$reference_genome, data = x$genomic_coordinates) %>% dplyr::filter(chr %in% chromosomes)

  low = min(reference_genome$from)
  upp = max(reference_genome$to)

  # Default blank genome -- remove labels with label_chr = NA
  bl_genome = suppressMessages(
    CNAqc:::blank_genome(label_chr = NA, ref = x$reference_genome, genomic_coords = x$genomic_coordinates) +
      labs(x = "", y = "")
  )

  # Segment id for the y-axis
  # bl_genome =
  bl_genome +
    ggplot2::geom_segment(
      data = calls_flat,
      ggplot2::aes(
        x = from,
        xend = to,
        color = label
      ),
      y = 1,
      yend = 1,
      size = 10
    ) +
    ggplot2::scale_color_manual(values = KARYO_colors) +
    ggplot2::coord_polar(theta = 'x', clip = 'off') +
    ggplot2::guides(color = F) +
    # guides(color = guide_legend('Karyotype', nrow = 1)) +
    ggplot2::ylim(-2, 2) +
    ggplot2::theme_void()
}
