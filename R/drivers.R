has_driver_data = function(x)
{
  stopifnot(inherits(x, 'cnaqc'))

  cn = colnames(x$mutations)

  if (all(c("is_driver", "driver_label") %in% cn)) return(TRUE)

  if ("is_driver" %in% cn & !("driver_label" %in% cn)){
    cli::cli_warn("Column 'is_driver' is annotated but 'driver_label' no -- did you try to add drivers data?")
  }

  if ("driver_label" %in% cn & !("is_driver" %in% cn)){
    cli::cli_warn("Column 'driver_label' is annotated but 'is_driver' no -- did you try to add drivers data?")
  }

  return(FALSE)
}

has_drivers = function(x)
{
  if(!has_driver_data(x)) return(FALSE)

  x$mutations %>%
    dplyr::filter(is_driver) %>%
    nrow() > 0
}

#' Extract drivers data.
#'
#' @description
#'
#' Returns a tibble with the drivers annotated on certain chromosomes,
#' extracting VAF or CCF values (if computed). If drivers are not
#' annotated \code{NULL} is returned.
#'
#' @param x A CNAqc object.
#' @param chromosomes Which chromosomes to subset.
#' @param which A keyword for \code{"VAF"} or \code{"CCF"}.
#'
#' @return A tibble
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(mutations = example_dataset_CNAqc$mutations, cna = example_dataset_CNAqc$cna, purity = example_dataset_CNAqc$purity)
#'
#' get_drivers(x)
get_drivers = function(x,
                       chromosomes = paste0('chr', c(1:22, 'X', 'Y')),
                       which = 'VAF')
{
  if(!has_driver_data(x)) return(NULL)

  # CCF?
  if (all(is.null(x$CCF_estimates)) & which == "CCF")
  {
    warning("Input does not have CCF estimates, see ?compute_CCF to determine CCF values.")
    return(NULL)
  }

  drivers_list = NULL

  if(which == "VAF")
    drivers_list = x$mutations %>%
      dplyr::filter(is_driver, chr %in% chromosomes)

  if(which == "CCF")
    drivers_list = CNAqc::CCF(x) %>%
      dplyr::filter(is_driver, chr %in% chromosomes)

  return(drivers_list)
}

# Add driver annotation to a histogram
annotate_drivers_to_histogram = function(x, drivers_list, p, which)
{
  # Annotate driver events if they exist
  if(is.null(drivers_list)) return(p)
  if(nrow(drivers_list) == 0) return(p)

  # Coordinate of the plot, place label in top part
  L = ggplot2::ggplot_build(p)$layout$panel_params[[1]]

  drivers_list$y = L$y.range[2] * .7

  # This overwrites the fill in p, is generally wrong
  # p =
  #   p +
  #   geom_vline(
  #     data = drivers_list,
  #     show.legend = FALSE,
  #     aes(color = karyotype, xintercept = eval(parse(text = which))),
  #     linetype = 'dashed',
  #     size = .3
  #   ) +
  #   scale_color_manual(values = CNAqc:::get_karyotypes_colors(unique(drivers_list$karyotype))) +
  #   geom_text(
  #     data = drivers_list,
  #     aes(
  #       x = eval(parse(text = which)),
  #       y = y,
  #       label = driver_label,
  #       fill = karyotype
  #     ),
  #     ylim = c(L$y.range[2] * .9, NA),
  #     size = 2,
  #     nudge_y = 0,
  #     nudge_x = 0,
  #     show.legend = FALSE
  #   ) +
  #     scale_fill_manual(values = CNAqc:::get_karyotypes_colors(unique(drivers_list$karyotype))) +
  #     coord_cartesian(clip = 'off')

  p = p +
    ggplot2::geom_vline(
      data = drivers_list,
      show.legend = FALSE,
      ggplot2::aes(color = karyotype, xintercept = eval(parse(text = which))),
      linetype = 'dashed',
      size = .3
    ) +
    ggplot2::scale_color_manual(values = CNAqc:::get_karyotypes_colors(unique(drivers_list$karyotype))) +
    ggrepel::geom_label_repel(
      data = drivers_list,
      ggplot2::aes(
        x = eval(parse(text = which)),
        y =  L$y.range[2] * .5,
        label = driver_label,
      ),
      color = 'black',
      # ylim = c(0.7, L$y.range[2] * .9),
      ylim = c(L$y.range[2] * .5, L$y.range[2]),
      size = 2,
      nudge_x = 0,
      show.legend = FALSE
    ) +
  ggplot2::coord_cartesian(clip = 'off') +
  ggplot2::labs(y = 'count')



  return(p)
}

# Add driver annotation to a segments plot
add_drivers_to_segment_plot = function(x, drivers_list, base_plot)
{
  # Annotate driver events if required
  if(is.null(drivers_list)) return(base_plot)
  if (nrow(drivers_list) == 0) return(base_plot)

  # Coordinate of the plot, place label in top part
  L = ggplot2::ggplot_build(base_plot)$layout$panel_params[[1]]

  drivers_list = CNAqc:::relative_to_absolute_coordinates(
    x,
    drivers_list %>% dplyr::filter(is_driver)
  )

  drivers_list$y = L$y.range[2] * .9

  base_plot +
    ggplot2::geom_vline(
      data = drivers_list,
      show.legend = FALSE,
      ggplot2::aes(xintercept = from),
      linetype = 'dashed',
      color = 'black',
      size = .3
    ) +
    ggrepel::geom_label_repel(
      data = drivers_list,
      ggplot2::aes(
        x = from,
        y = y,
        label = driver_label,
        fill = karyotype
      ),
      ylim = c(L$y.range[2] * .9, NA),
      size = 2,
      nudge_y = 0,
      nudge_x = 0,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = CNAqc:::get_karyotypes_colors(unique(drivers_list$karyotype))) +
    ggplot2::coord_cartesian(clip = 'off')



}
