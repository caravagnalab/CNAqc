#' Title
#'
#' @param x
#' @param chromosomes
#'
#' @return
#' @export
#'
#' @examples
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



