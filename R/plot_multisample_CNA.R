# data('example_evoverse', package = 'evoverse')
# x = example_evoverse
# L = x$CNAqc



#' Plots CNA for multiple samples.
#'
#' @description
#'
#' This functions plots with circular layout the CNA segments of multiple
#' `CNAqc` objects, stored in a named list `x`. In this layout, instead of
#' plotting the Major and minor alleles of each segment, we assign colours
#' to certain karyotypes. Each sample is plot on a lane, like in a donut plot.
#'
#' @param x A named list of `CNAqc` objects.
#'
#' @return A `ggplot` plot
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna,example_dataset_CNAqc$purity)
#'
#' plot_multisample_CNA(list(`S1` = x, `S2` = x))
plot_multisample_CNA = function(x)
{
  L = x
  Ln = names(L)
  if(is.null(Ln)) Ln = paste0("Sample ", 1:length(L))

  KARYO_colors =
    c(
      `0:0` = 'darkblue',
      `1:0` = 'steelblue',
      `1:1` = ggplot2::alpha('forestgreen', .8),
      `2:0` = 'turquoise4',
      `2:1` = ggplot2::alpha('orange', .8),
      `2:2` = 'indianred3',
      `other` = 'darkred'
    )

  # Extract calls, and flatten them for plotting
  calls = lapply(Ln,
                 function(s)
                 {
                   L[[s]]$cna %>%
                     mutate(
                       label = paste(Major, minor, sep = ':'),
                       CN = minor + Major,
                       sample = s
                     ) %>%
                     select(chr, from, to, label, CN, sample)
                 })

  calls = lapply(calls, CNAqc:::relative_to_absolute_coordinates)

  calls_flat =
    suppressWarnings(
      Reduce(
        function(x, y) full_join(x, y, by = c("chr", "from", "to", "label", "CN", "sample")),
        calls) %>%
      mutate(
        label = ifelse(label %in% names(KARYO_colors), label, 'other')
      )
    )

  chromosomes = calls_flat$chr %>% unique

  # Get hg19 coordinates for used chromosomes
  data('chr_coordinates_hg19', package = 'CNAqc')

  chr_coordinates_hg19 = chr_coordinates_hg19 %>% filter(chr %in% chromosomes)

  low = min(chr_coordinates_hg19$from)
  upp = max(chr_coordinates_hg19$to)

  # Default blank genome -- remove labels with label_chr = NA
  bl_genome = suppressMessages(
    CNAqc:::blank_genome(label_chr = NA) +
      labs(x = "", y = "")
  )

  # Segment id for the y-axis
  seg_id = pio:::nmfy(Ln, seq_along(Ln))
  calls_flat$sample_id = seg_id[calls_flat$sample]

  # bl_genome =
  bl_genome +
    geom_segment(
      data = calls_flat,
      aes(
        x = from,
        xend = to,
        y = sample_id,
        yend = sample_id,
        color = label
      ),
      size = 5
    ) +
    scale_color_manual(values = KARYO_colors) +
    coord_polar(theta = 'x', clip = 'off') +
    guides(color = guide_legend('Karyotype', nrow = 1)) +
    ylim(-5, max(seg_id) + 3) +
    labs(
      title = "Comparative CNA",
      subtitle = paste0('Tracks: ', paste(Ln, collapse = ', '))
    ) +
    theme(
      legend.key.height = unit(.1, "cm"),
      axis.text.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(size = .3)
    )

}


