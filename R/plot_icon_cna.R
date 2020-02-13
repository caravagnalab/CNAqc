#' Plot icon for these segments.
#'
#' @description
#'
#' An icon is a small plot with no legend or axis, which
#' has the same layout of function \code{plot_multisample_CNA}.
#'
#' @param x A `CNAqc` object.
#'
#' @return A `ggplot` plot.
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna,example_dataset_CNAqc$purity)
#'
#' plot_icon_CNA(x)
#'
#' # Same legend of this
#' plot_multisample_CNA(list(`S1` = x, `S2` = x))
plot_icon_CNA = function(x)
{
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
  reference_genome = CNAqc:::get_reference(x$reference_genome) %>% dplyr::filter(chr %in% chromosomes)

  low = min(reference_genome$from)
  upp = max(reference_genome$to)

  # Default blank genome -- remove labels with label_chr = NA
  bl_genome = suppressMessages(
    CNAqc:::blank_genome(label_chr = NA) +
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
