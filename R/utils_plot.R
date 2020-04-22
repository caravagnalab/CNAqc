my_ggplot_theme = function(cex = 1)
{
  cex_opt = getOption('CNAqc_cex', default = 1)

  theme_light(base_size = 10 * cex_opt) +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(.3 * cex_opt, "cm"),
      panel.background = element_rect(fill = 'white')
    )
}

get_karyotypes_colors = function(karyotypes)
{
  color = c(
    `1:1` = ggplot2::alpha('forestgreen', .8),
    '1:0' = 'steelblue',
    '0:0' = 'darkblue',
    '2:0' = 'turquoise4',
    '2:1' = ggplot2::alpha('orange', .8),
    '2:2' = 'firebrick3'
  )

  missing = setdiff(karyotypes, names(color))
  nmissing = length(missing)


  c(color, pio:::nmfy(missing, rep('gray', nmissing)))
}

add_segments_to_plot = function(segments, base_plot, cn)
{
  if (cn == 'absolute')
  {
    # Add one Major and minor lines, one on top of the other, red and blu
    M_seg = segments %>% dplyr::select(from, to, Major) %>% dplyr::rename(value = Major)
    m_seg = segments %>% dplyr::select(from, to, minor) %>% dplyr::rename(value = minor)

    base_plot = base_plot +
      geom_segment(
        data = M_seg %>% dplyr::mutate(Allele = "Major allele"),
        aes(
          x = from,
          xend = to,
          y = value,
          yend = value,
          colour = Allele
        ),
        size = 1.5
      ) +
      geom_segment(
        data = m_seg %>% dplyr::mutate(Allele = "minor allele"),
        aes(
          x = from,
          xend = to,
          y = value - 0.1,
          yend = value - 0.1,
          colour = Allele
        ),
        size = 1
      ) +
      scale_color_manual(values = c(`Major allele` = 'red', `minor allele` = 'steelblue')) +
      guides(color = guide_legend(''))

    # Some layout
    base_plot = base_plot +
      theme(
        legend.position = "bottom",
        legend.justification = "right",
        legend.margin = margin(0, 0, 0, 0)
      ) +
      labs(y = "Absolute allele counts")

  }

  if (cn == 'total')
  {
    t_seg = segments %>% dplyr::select(from, to, total) %>% rename(value = total)

    base_plot = base_plot +
      geom_segment(
        data = t_seg %>% dplyr::mutate(Allele = "Segment ploidy"),
        aes(
          x = from,
          xend = to,
          y = total,
          yend = total
        ),
        size = 1.5,
        colour = 'black'
      )
  }

  return(base_plot)
}

add_shadow_to_plot = function(segments, base_plot,  highlight)
{
  # Shadow CN segments, if any
  all_karyotypes = CNAqc:::get_karyotypes_colors(NULL)
  all_karyotypes = all_karyotypes[names(all_karyotypes) %in% highlight]

  segments = segments %>%
     dplyr::filter(karyotype %in% names(all_karyotypes))

  if (nrow(segments) > 0)
    base_plot = base_plot +
    geom_rect(
      data = segments,
      aes(
        xmin = from,
        xmax = to,
        ymin = -Inf,
        ymax = Inf,
        fill = factor(karyotype, levels = c('2:0', '1:0', '1:1', '2:1', '2:2'))
      ),
      alpha = .3
    ) +
    scale_fill_manual(values = all_karyotypes) +
    guides(fill = guide_legend('', override.aes = list(alpha = 1)))

  return(base_plot)
}

add_breakpoints_to_plot = function(segments, base_plot, max_Y_height)
{
  # Capped off segments too high to render the plot readable
  off_plot = segments %>% dplyr::filter(total > max_Y_height)

  # If any, add a top-label and a line
  if (nrow(off_plot) > 0)
  {
    base_plot = base_plot +
      geom_hline(
        yintercept = max_Y_height,
        size = .2,
        color = 'darkgray',
        linetype = 'dashed'
      )

    # Simulate an internal legendq
    L = ggplot_build(base_plot)$layout$panel_params[[1]]
    Lx = abs(L$x.range[2] - L$x.range[1]) * .85

    base_plot = base_plot +
      geom_label(
        data = data.frame(
          x = Lx,
          y = L$y.range[2] - 0.5,
          label = paste0('< ', max_Y_height)
        ),
        aes(x = x, y = y, label = label),
        fill = 'darkgray',
        color = 'white',
        size = 2,
        nudge_y = 0,
        nudge_x = 0,
        inherit.aes = FALSE
      )
  }

  # Minimum height of the plot
  L = ggplot_build(base_plot)$layout$panel_params[[1]]

  if (L$y.range[2] <= 5) {
    base_plot = base_plot + ylim(-0.5, 5)
  }


  # =-=-=-=-=-=-=-=-=-=-=-=-
  # Breakpoints annotations
  # =-=-=-=-=-=-=-=-=-=-=-=-
  breakpoints = data.frame(
    x = segments$from,
    y = 0.1,
    outern = segments$Major > max_Y_height
  )

  base_plot = base_plot +
    geom_point(
      data = breakpoints %>% filter(!outern),
      aes(x = x, y = y),
      size = .5,
      shape = 1,
      color = 'darkgray'
    ) +
    geom_point(
      data = breakpoints %>% filter(outern),
      aes(x = x, y = y),
      size = .5,
      color = 'black'
    )

  return(base_plot)
}

add_drivers_to_segment_plot = function(x, drivers_list, base_plot)
{
  # Annotate driver events if required
  if (!("is_driver" %in% colnames(x$snvs))) return(base_plot)

  drivers_list = drivers_list %>% dplyr::filter(is_driver)
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
        # fill = karyotype
      ),
      ylim = c(L$y.range[2] * .9, NA),
      size = 2,
      nudge_y = 0,
      nudge_x = 0,
      show.legend = FALSE
    )


}

relative_to_absolute_coordinates = function(x, cna)
{
  reference_genome = CNAqc:::get_reference(x$reference_genome)

  vfrom = reference_genome$from
  names(vfrom) = reference_genome$chr

  cna %>%
    mutate(from = from + vfrom[chr],
           to = to + vfrom[chr])
}

absolute_to_relative_coordinates = function(x, cna)
{
  reference_genome = CNAqc:::get_reference(x$reference_genome)

  vfrom = reference_genome$from
  names(vfrom) = reference_genome$chr

  cna %>%
    mutate(from = from - vfrom[chr],
           to = to - vfrom[chr])
}
