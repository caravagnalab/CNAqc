#' Plot CNA segments.
#'
#' @description Plot CNA segments as Major/ minor allele, annotating
#' clonal and subclonal CNA calls in two different sets of colors.
#' This function uses \code{hg19} genome coordinates to map CNA segments.
#'
#' @param x An object of class \code{cnaqc}, created by the \code{init} function.
#' @param chromosomes The chromosome to use for this plot.
#' @param max_Y_height Maximum height onn the Y-axis of the plot, if there are segments
#' above this heigh a warning is raised and annnotated in the plot as well in graphical format.
#' @param circular Uses a circular layout in polar coordinates to make the segments
#' look like a circos plot. It can save space.
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
                         max_Y_height = 8,
                         circular = FALSE,
                         ...)
{
  stopifnot(inherits(x, 'cnaqc'))

  # Circular layout
  if(circular) return(plot_segments_circular(x, chromosomes = chromosomes))

  # Standard plot -- baseline genome reference
  base_plot = CNAqc:::blank_genome(chromosomes = chromosomes, ref = x$reference_genome, ...)

  # Segments
  segments = x$cna %>%
    filter(chr %in% chromosomes)
  segments = CNAqc:::relative_to_absolute_coordinates(x, segments)

  # if there are 1-1 segments, shadow them
  one_one = segments %>% filter(Major == 1, minor == 1)

  if (nrow(one_one) > 0)
    base_plot = base_plot +
    geom_rect(
      data = one_one,
      aes(
        xmin = from,
        xmax = to,
        ymin = -Inf,
        ymax = Inf
      ),
      alpha = .2,
      fill = 'forestgreen'
    )

  # Fragmentation ~ add some annotation to hihglight that
  if(!is.null(x$arm_fragmentation))
  {
    fragmented = x$arm_fragmentation$table %>%
      dplyr::filter(significant, chr %in% chromosomes) %>%
      dplyr::mutate(label = paste0(chr, arm)) %>%
      dplyr::pull(label)

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

  # Segments
  base_plot = base_plot +
    geom_segment(
      data = segments,
      aes(
        x = from,
        xend = to,
        y = Major,
        yend = Major
      ),
      size = 1.5,
      colour = 'red'
    ) +
    geom_segment(
      data = segments %>% mutate(minor = minor - 0.1),
      aes(
        x = from,
        xend = to,
        y = minor,
        yend = minor
      ),
      size = 1,
      colour = 'steelblue'
    ) +
    labs(
      caption =
        paste0(
          x$n_snvs, ' mutations, ',
          x$n_cna_clonal, ' clonal CNA, ',
          x$n_cna_sbclonal, ' subclonal CNA'
        )
    )

  # =-=-=-=-=-=-=-=-=-=-=-=-
  # Be smart layout
  # =-=-=-=-=-=-=-=-=-=-=-=-
  # 1) chop off segments too high to render the plot readable
  MH = max(max(segments$minor), max(segments$Major))

  if(MH > max_Y_height) {
    warning("Segments with CN above ", max_Y_height, " will not be plot; this is annotated in the figure.")

    base_plot = base_plot + ylim(-0.5, max_Y_height)
  }

  # 2) minimum height of the plot
  if(MH <= max_Y_height) {
    warning("No segments with CN above 3, the plot is anyway scaled up to ", max_Y_height, " to render better.")

    base_plot = base_plot + ylim(-0.5, 5)
  }

  # =-=-=-=-=-=-=-=-=-=-=-=-
  # Breakpoints annotations
  # =-=-=-=-=-=-=-=-=-=-=-=-
  breakpoints = data.frame(x = segments$from, y = 0.1, outern = segments$Major > max_Y_height)

  base_plot = base_plot +
    geom_point(
      data = breakpoints %>% filter(!outern),
      aes(x = x, y = y),
      size = .5,
      shape = 1,
      color = 'black'
    ) +
    geom_point(
      data = breakpoints %>% filter(outern),
      aes(x = x, y = y),
      size = .5,
      shape = 4,
      color = 'orange'
    )

  # Annotate driver events if required
  if("is_driver" %in% colnames(x$snvs))
  {
    driver_list = x$snvs %>%
      dplyr::filter(is_driver, chr %in% chromosomes) %>%
      dplyr::left_join(
        x$cna %>% dplyr::select(segment_id, Major),
        by = 'segment_id')

    if(nrow(driver_list) > 0)
      base_plot = base_plot +
        ggrepel::geom_text_repel(
          data = CNAqc:::relative_to_absolute_coordinates(x, driver_list),
          aes(x = from, y = Major, label = gene),
          ylim = c(max_Y_height - 1, NA),
          segment.colour = 'darkblue',
          segment.size = .3,
          nudge_y = max_Y_height - 0.5,
          color = 'darkblue',
          arrow = grid::arrow(length = unit(0.05, "inches"), type  = 'closed')
        )
  }

  return(base_plot)
}


# Internal function -- implements a circular layout using plot_segments
plot_segments_circular = function(x, chromosomes = paste0('chr', c(1:22, 'X', 'Y')))
{
  # Uses the standard function, but rotates the coordinate system,
  # offsets the y-axis and add a new set of labels
 suppressMessages(
   plot_segments(x, max_Y_height = 5, chromosomes = chromosomes, circular = FALSE, label_chr = NA) +
      coord_polar(theta = 'x', start = 0, clip = 'off') +
      ylim(-2, 5) +
     labs(
       x = "",
       y = ""
     ) +
     theme(
       axis.text.y = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       panel.border = element_rect(size = .3)
     )
 )
}


