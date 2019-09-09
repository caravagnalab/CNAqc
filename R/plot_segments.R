#' Plot CNA segments.
#'
#' @description Plot CNA segments as Major/ minor allele, annotating
#' clonal and subclonal CNA calls in two different sets of colors.
#' This function uses \code{hg19} genome coordinates to map CNA segments.
#'
#' @param x An object of class \code{cnaqc}, created by the \code{init} function.
#' @param chromosomes The chromosome to use for this plot.
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
plot_segments = function(x, chromosomes = paste0('chr', c(1:22, 'X', 'Y')))
{

  stopifnot(inherits(x, 'cnaqc'))

  base_plot = blank_genome(chromosomes)

  # Segments
  segments = x$cna %>%
    filter(chr %in% chromosomes) %>%
    relative_to_absolute_coordinates

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

  base_plot
}



