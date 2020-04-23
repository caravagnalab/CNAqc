has_driver_data = function(x)
{
  if (all(c("is_driver", "driver_label") %in% colnames(x$snvs))) return(TRUE)

  return(FALSE)
}

#' Extract driver data.
#'
#' @description
#'
#' Returns a tibble with the drivers annotated on certain chromsomes,
#' extracting VAF or CCF values (if computed). If drivers are not
#' annotated \code{NULL} is returned.
#'
#' @param x An object of class \code{cnaqc}, created by the \code{init} function
#' that can have drivers annotated.
#' @param chromosomes Which chromosome to subset
#' @param which A keyword for \code{"VAF"} or \code{"CCF"}.
#'
#' @return A tibble
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna,example_dataset_CNAqc$purity)
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
    drivers_list = x$snvs %>%
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
  L = ggplot_build(p)$layout$panel_params[[1]]

  drivers_list$y = L$y.range[2] * .9

  p = p +
    geom_vline(
      data = drivers_list,
      show.legend = FALSE,
      aes(color = karyotype, xintercept = eval(parse(text = which))),
      linetype = 'dashed',
      size = .3
    ) +
    scale_color_manual(values = CNAqc:::get_karyotypes_colors(unique(drivers_list$karyotype))) +
    ggrepel::geom_label_repel(
      data = drivers_list,
      aes(
        x = eval(parse(text = which)),
        y = y,
        label = gene,
        fill = karyotype
      ),
      ylim = c(L$y.range[2] * .9, NA),
      size = 2,
      nudge_y = 0,
      nudge_x = 0,
      show.legend = FALSE
    )

  return(p)
}

# Add driver annotation to a segments plot
add_drivers_to_segment_plot = function(x, drivers_list, base_plot)
{
  # Annotate driver events if required
  if (nrow(drivers_list) == 0) return(base_plot)

  # Coordinate of the plot, place label in top part
  L = ggplot_build(base_plot)$layout$panel_params[[1]]

  drivers_list = CNAqc:::relative_to_absolute_coordinates(
    x,
    drivers_list %>% dplyr::filter(is_driver)
  )

  drivers_list$y = L$y.range[2] * .9

  base_plot +
    geom_vline(
      data = drivers_list,
      show.legend = FALSE,
      aes(xintercept = from),
      linetype = 'dashed',
      color = 'black',
      size = .3
    ) +
    ggrepel::geom_label_repel(
      data = drivers_list,
      aes(
        x = from,
        y = y,
        label = gene,
        fill = karyotype
      ),
      ylim = c(L$y.range[2] * .9, NA),
      size = 2,
      nudge_y = 0,
      nudge_x = 0,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = CNAqc:::get_karyotypes_colors(NULL))


}
