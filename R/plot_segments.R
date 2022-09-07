#' Plot CNA segments.
#'
#' @description
#' CNAqc can plot genome-wide segments showing major and minor allele counts as
#' red and blue bars. Bottom circles annotate breakpoints; by default this plot
#' has limited y-axis and breakpoints for segments outside the plot (e.g. very
#' high amplifications) are in black. Areas in the genome that are mapped to
#' the most prevalent karyotype are shadowed by default. Note that the colour
#' scheme in CNAqc is fixed for certain segments. The colour scheme is adopted
#' also for other information, e.g., to report the segment where a certain driver
#' is mapped. A circular layout for this plot is also available, and total
#' copy numbers can also be plot.
#'
#' Note that two distinct colour schemas distinguish clonal and subclonal segments.
#'
#' @param x A CNAqc object.
#' @param chromosomes The chromosome to plot.
#' @param max_Y_height Maximum height for the Y-axis of the plot. Segments witht total copy
#' number (Major plus minor) above \code{max_Y_height} are not showm; in that case a point
#' annotation is put on top of the plot, and \code{max_Y_height} is annotated in the grid.
#' @param highlight A list of karyotype ids in \code{"Major:minor"} notation
#' (e.g., \code{"1:1", "2,1", ...}) that will be shadowed with a transparent area.
#' By default, it plots the most prevalent karyotype using `x$most_prevalent_karyotype`.
#' @param circular Uses a circular layout in polar coordinates to make the segments
#' look like a circos plot. This visualisation can save space.
#' @param cn Type of copy number segment to show on the plot. Either \code{"absolute"} for
#' Major and minor annotation, or \code{"total"} for the total (Major plus minor) annotation.
#'
#' @return A \code{ggplot2} plot.
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(mutations = example_dataset_CNAqc$mutations, cna = example_dataset_CNAqc$cna, purity = example_dataset_CNAqc$purity)
#'
#' # Default view, works for most cases
#' plot_segments(x)
#'
#' # A subset of chromosomes
#' plot_segments(x, chromosomes = 'chr13')
#'
#' # Total CNAs
#' plot_segments(x, cn = 'total')
#'
#' # Shadows all simple clonal CNAs
#' plot_segments(x, highlight = c('1:0', '1:1', '2:0', '2:1', '2:2'))
plot_segments = function(x,
                         chromosomes = paste0('chr', c(1:22, 'X', 'Y')),
                         max_Y_height = 6,
                         circular = FALSE,
                         cn = 'absolute',
                         highlight = x$most_prevalent_karyotype,
                         ...)
{
  stopifnot(inherits(x, 'cnaqc'))

  # Circular layout
  if (circular) {
    base_plot = plot_segments_circular(x, chromosomes = chromosomes)

    return(base_plot)
  }
  else{
    # Standard plot -- baseline genome reference
    base_plot = blank_genome(chromosomes = chromosomes,
                                     ref = x$reference_genome,
                                     ...)
  }

  # Segments
  segments = x$cna %>%
    dplyr::filter(chr %in% chromosomes)  %>%
    dplyr::mutate(
      total = Major + minor,
      karyotype = paste0(Major, ':', minor)
      )

  segments = CNAqc:::relative_to_absolute_coordinates(x, segments)

  # =-=-=-=-=-=-=-=-=-=-=-=-
  # Shadow for highligthing
  # =-=-=-=-=-=-=-=-=-=-=-=-
  base_plot = CNAqc:::add_shadow_to_plot(segments, base_plot, highlight)

  # =-=-=-=-=-=-=-=-=-=-=-=-
  # Draw Segments
  # =-=-=-=-=-=-=-=-=-=-=-=-
  base_plot = add_segments_to_plot(
    segments = segments %>% dplyr::filter(total <= max_Y_height),
    base_plot = base_plot,
    cn = cn)

  # Extract subclonal segments
  subclonal_segments = NULL
  if (!is.null(x$cna_subclonal) & nrow(x$cna_subclonal) > 0)
  {
    subclonal_segments = x$cna_subclonal %>%
      dplyr::filter(chr %in% chromosomes)

    if (nrow(subclonal_segments) > 0)
    {
      base_plot = add_subclonal_segments_to_plot(
        segments = subclonal_segments %>%
          relative_to_absolute_coordinates(x = x),
        base_plot = base_plot,
        cn = cn
      )
    }
  }

  # Fragmentation ~ add some annotation to hihglight that
  if (!is.null(x$arm_fragmentation))
  {
    fragmented = x$arm_fragmentation$table %>%
      dplyr::filter(significant, chr %in% chromosomes) %>%
      dplyr::mutate(label = paste0(chr, arm)) %>%
      dplyr::pull(label)

    if (length(fragmented) > 0)
    {
      expanded_reference = CNAqc:::expand_reference_chr_to_arms(x) %>%
        dplyr::filter(chr %in% fragmented)

      # base_plot = base_plot +
      #   geom_rect(
      #     data = expanded_reference,
      #     aes(
      #       xmin = from,
      #       xmax = to,
      #       ymin = -Inf,
      #       ymax = Inf
      #     ),
      #     alpha = .2,
      #     fill = NA,
      #     color = 'purple4'
      #   )
      #
      # base_plot +
      #   geom_segment(
      #     data = expanded_reference,
      #     aes(
      #       x = from,
      #       xend = to,
      #       y = -0.2,
      #       yend = -0.2,
      #     ),
      #     color = 'purple4',
      #     linetype = 1,
      #     size = 2
      #   )

      base_plot = base_plot +
        ggplot2::geom_label(
          data = expanded_reference,
          ggplot2::aes(
            x = from,
            label = gsub(pattern = 'chr', replacement = '', chr),
            y = -0.2
          ),
          fill = 'purple4',
          color = 'white',
          size = 2,
          label.padding = unit(0.05, 'cm')
        )
      }
  }

  # Caption

  capt_label = paste0(
    "Ploidy ", x$ploidy, "; Purity  ", x$purity,
    '; n = ', x$n_mutations, ' mutations in ',
    x$n_cna_clonal,
    ' segments'
  )

  # Add extras
  if (!is.null(x$arm_fragmentation))
    capt_label = paste0(
      capt_label,
      '; ',
      x$arm_fragmentation$table %>%
        dplyr::filter(significant, chr %in% chromosomes) %>%
        length,
      ' fragmented arms'
    )

  base_plot = base_plot + ggplot2::labs(caption = capt_label)

  # =-=-=-=-=-=-=-=-=-=-=-=-
  # Breakpoints annotations
  # =-=-=-=-=-=-=-=-=-=-=-=-
  base_plot = CNAqc:::add_breakpoints_to_plot(segments, base_plot, max_Y_height, circular)

  # =-=-=-=-=-=-=-=-=-=-=-=-
  # Drivers annotations
  # =-=-=-=-=-=-=-=-=-=-=-=-
  drivers_list = CNAqc:::get_drivers(x, chromosomes = chromosomes)
  if(!circular)
  base_plot = CNAqc:::add_drivers_to_segment_plot(x, drivers_list = drivers_list, base_plot)


  return(base_plot)
}


# Internal function -- implements a circular layout using plot_segments
plot_segments_circular = function(x, chromosomes = paste0('chr', c(1:22, 'X', 'Y')))
{
  # Uses the standard function, but rotates the coordinate system,
  # offsets the y-axis and add a new set of labels
  suppressMessages(
    plot_segments(
      x,
      max_Y_height = 5,
      chromosomes = chromosomes,
      circular = FALSE,
      label_chr = NA
    ) +
      ggplot2::coord_polar(
        theta = 'x',
        start = 0,
        clip = 'off'
      ) +
      ggplot2::ylim(-2, 5) +
      ggplot2::labs(x = "",
           y = "") +
      ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(size = .3)
      )
  )
}
