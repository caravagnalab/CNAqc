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

  missing = setdiff(names(color),karyotypes)
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
        data = M_seg %>% dplyr::mutate(Allele = "Major allele (clonal)"),
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
        data = m_seg %>% dplyr::mutate(Allele = "minor allele (clonal)"),
        aes(
          x = from,
          xend = to,
          y = value - 0.1,
          yend = value - 0.1,
          colour = Allele
        ),
        size = 1
      ) +
      scale_color_manual(values = c(`Major allele (clonal)` = 'red', `minor allele (clonal)` = 'steelblue')) +
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
    base_plot = base_plot +
      geom_segment(
        data = segments %>% dplyr::select(from, to, total) %>% dplyr::mutate(Allele = "Segment ploidy"),
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

add_subclonal_segments_to_plot = function(segments, base_plot, cn)
{
  if (cn == 'absolute')
  {
    # Add one Major and minor lines per clone
    M_seg_1 = segments %>% dplyr::select(from, to, Major, CCF) %>% dplyr::rename(value = Major)
    m_seg_1 = segments %>% dplyr::select(from, to, minor, CCF) %>% dplyr::rename(value = minor)
    M_seg_2 = segments %>% dplyr::select(from, to, Major_2, CCF) %>% dplyr::rename(value = Major_2) %>% mutate(CCF = 1 - CCF)
    m_seg_2 = segments %>% dplyr::select(from, to, minor_2, CCF) %>% dplyr::rename(value = minor_2) %>% mutate(CCF = 1 - CCF)

    sbc = RColorBrewer::brewer.pal(10, "Paired")
    sbc = sbc[7:10]

    # Re-define colours also for clonal ones
    my_colors = c(
      `Major allele (clonal)` = 'red',
      `minor allele (clonal)` = 'steelblue',
      `Major allele (Subclone 1)` = sbc[2],
      `minor allele (Subclone 1)` = sbc[1],
      `Major allele (Subclone 2)` = sbc[4],
      `minor allele (Subclone 2)` = sbc[3]
    )

    # linetypes = c(
    #   `<33%` = "CCF <= 33%",
    #   `33%-66%` = "33% < CCF <= 66%",
    #   `>66%` = "CCF > 66%"
    # )

    # Map colours
    base_plot = suppressWarnings(suppressMessages(
      base_plot +
        scale_color_manual(values = my_colors)+
        guides(linetype = guide_legend("Subclone CCF"))
    ))

    for(i in 1:nrow(segments))
    {
      M1 = M_seg_1 %>% dplyr::mutate(Allele = "Major allele (Subclone 1)") %>% filter(row_number() == i)
      m1 = m_seg_1 %>% dplyr::mutate(Allele = "minor allele (Subclone 1)") %>% filter(row_number() == i)
      M2 = M_seg_2 %>% dplyr::mutate(Allele = "Major allele (Subclone 2)") %>% filter(row_number() == i)
      m2 = m_seg_2 %>% dplyr::mutate(Allele = "minor allele (Subclone 2)") %>% filter(row_number() == i)

      # Shift height to avoid overlap
      i_segments = bind_rows(M1, m1, M2, m2) %>%
        group_by(value) %>%
        mutate(pos = row_number() - 1) %>%
        mutate(value = value - 0.1*pos)
      # %>%
      #   mutate(linetype = case_when(
      #     CCF <= .33 ~ linetypes[1],
      #     CCF > .33 & CCF <= .66 ~ linetypes[2],
      #     CCF > .66 ~ linetypes[3]
      #   ))

      base_plot = base_plot +
        geom_segment(
          data = i_segments %>% filter(grepl("Major", Allele)),
          aes(
            x = from,
            xend = to,
            y = value,
            yend = value,
            colour = Allele
            # linetype = linetype %>% paste()
          ),
          size = 1.5
        ) +
        geom_segment(
          data = i_segments %>% filter(grepl("minor", Allele)),
          aes(
            x = from,
            xend = to,
            y = value,
            yend = value,
            colour = Allele
            # linetype = linetype %>% paste()
          ),
          size = 1
        )
    }

  }

  if (cn == 'total')
  {
    message("TODO - yet to implement.")
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

add_breakpoints_to_plot = function(segments, base_plot, max_Y_height, circular)
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

  if(!circular){
    if (L$y.range[2] <= 5) {
      base_plot = base_plot + ylim(-0.5, 5)
    }
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
