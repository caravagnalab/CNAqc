#' Plot CNA segments.
#'
#' @description Plot CNA segments as Major/ minor allele, annotating
#' clonal and subclonal CNA calls in two different sets of colors.
#' This function uses \code{hg19} genome coordinates to map CNA segments.
#'
#' @param x An object of class \code{cnaqc}, created by the \code{init} function.
#' @param chromosomes The chromosome to use for this plot.
#' @param max_Y_height Maximum height for the Y-axis of the plot. Segments witht total copy
#' number (Major plus minor) above \code{max_Y_height} are not showm; in that case a point
#' annotation is put on top of the plot, and \code{max_Y_height} is annotated in the grid.
#' @param highlight A list of karyotype ids in \code{"Major:minor"} notation
#' (e.g., \code{"1:1", "2,1", ...}) that will be shadowed with a transparent area.
#' By default, it plots the most prevalent karyotype,
#' @param circular Uses a circular layout in polar coordinates to make the segments
#' look like a circos plot. This visualisation can save space.
#' @param cn Type of copy number segment to show on the plot. Either \code{"absolute"} for
#' Major and minor annotation, or \code{"total"} for the total (Major plus minor) annotation.
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna,example_dataset_CNAqc$purity)
#'
#' plot_segments(x)
#' plot_segments(x, chromosomes = 'chr13')
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
  if (circular)
    return(plot_segments_circular(x, chromosomes = chromosomes))

  # Standard plot -- baseline genome reference
  base_plot = CNAqc:::blank_genome(chromosomes = chromosomes,
                                   ref = x$reference_genome,
                                   ...)

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
  base_plot = CNAqc:::add_segments_to_plot(
    segments = segments %>% dplyr::filter(total <= max_Y_height),
    base_plot = base_plot,
    cn = cn)

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

      base_plot = base_plot +
        geom_rect(
          data = expanded_reference,
          aes(
            xmin = from,
            xmax = to,
            ymin = -Inf,
            ymax = Inf
          ),
          alpha = .2,
          fill = NA,
          color = 'purple4'
        )

    }
  }

  # Caption
  base_plot = base_plot +
    labs(
      caption =
        paste0(
          "Ploidy ", x$ploidy, '; ',
          x$n_snvs,
          ' mutations, ',
          x$n_cna_clonal,
          ' clonal CNA, ',
          x$n_cna_sbclonal,
          ' subclonal CNA'
        )
    )

  # =-=-=-=-=-=-=-=-=-=-=-=-
  # Breakpoints annotations
  # =-=-=-=-=-=-=-=-=-=-=-=-
  base_plot = CNAqc:::add_breakpoints_to_plot(segments, base_plot, max_Y_height)

  # =-=-=-=-=-=-=-=-=-=-=-=-
  # Drivers annotations
  # =-=-=-=-=-=-=-=-=-=-=-=-
  drivers_list =  x$snvs %>% dplyr::filter(chr %in% chromosomes)
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
      coord_polar(
        theta = 'x',
        start = 0,
        clip = 'off'
      ) +
      ylim(-2, 5) +
      labs(x = "",
           y = "") +
      theme(
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = .3)
      )
  )
}
