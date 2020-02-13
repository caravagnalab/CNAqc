my_ggplot_theme = function(cex = 1)
{
  theme_light(base_size = 10 * cex) +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(.3 * cex, "cm"),
      panel.background = element_rect(fill = 'white')
    )
}

get_karyotypes_colors = function(karyotypes)
{
  color = c(
    '1:1' = 'forestgreen',
    '1:0' = 'steelblue',
    '0:0' = 'cadetblue',
    '2:0' = 'deepskyblue3',
    '2:1' = 'darkorange3',
    '2:2' = 'firebrick3'
  )

  missing = setdiff(karyotypes, names(color))
  nmissing = length(missing)


  c(color, pio:::nmfy(missing, rep('gray', nmissing)))
}

relative_to_absolute_coordinates = function(x, cna) 
{
  reference_genome = CNAqc:::get_reference(x$reference_genome)

  vfrom = reference_genome$from
  names(vfrom) = reference_genome$chr

  cna %>%
    mutate(
      from = from + vfrom[chr],
      to = to + vfrom[chr]
    )
}

