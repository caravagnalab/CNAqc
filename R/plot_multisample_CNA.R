#' Plots CNAs from multiple samples.
#'
#' @description
#'
#' This functions plots with distinct layout CNA data associated to multiple
#' CNAqc objects.
#'
#' Flat layout: classical layout where the amount of deletions and gains are
#' reported with a certain discretized binning of the input tumour genome. Deletions
#' are anything with an LOH state; gains must have more than 3 copies.
#'
#' Circular layout: plot the Major and minor alleles of each segment, each sample
#' is plot on a lane, like in a donut plot.
#'
#' @param x A list of `CNAqc` objects; names are optional for the circular layout.
#'
#' @return A `ggplot` plot.
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(example_dataset_CNAqc$mutations, example_dataset_CNAqc$cna, example_dataset_CNAqc$purity)
#'
#' # Add some example deletion
#' x2 = x
#' x2$cna$Major[1:10] = 2
#' x2$cna$minor[1:10] = 0
#'
#' plot_multisample_CNA(list(`S1` = x, `S2` = x))
#'
#'
#' plot_multisample_CNA(list(`S1` = x, `S2` = x), layout = 'circular')
plot_multisample_CNA = function(x, layout = 'flat', ...)
{
  ok_input = sapply(x, function(x) inherits(x, what = 'cnaqc')) %>% all()

  if(!ok_input) stop("Input x must be a list of CNAqc objects!")

  if(layout == "flat") return(x %>% aux_plot_cohort_CNA(...))
  if(layout == "circular") return(x %>% aux_plot_cohort_CNA_circular(...))

  return(ggplot())
}

# Classical CNA copy number plot
aux_plot_cohort_CNA = function(x, delta = 1e5)
{
  chop_clonal_segments = function(x, delta = 1e5)
  {
    # first, break the chromosomes into chunks of resolution at fixed size
    easypar::run(
      FUN = function(chr) {
        # Segments on chromomosme chr
        x_chr = x$cna %>% dplyr::filter(chr == !!chr)

        # Bin index
        x_chr = x_chr %>% dplyr::mutate(from_bin = floor(from / delta),
                                        to_bin = floor(to / delta))

        # Classify all bins
        x_chr = x_chr %>% dplyr::mutate(
          CNA = dplyr::case_when(
            minor == 0 & Major == 1 ~ "Deletion",
            minor == 0 & Major == 2 ~ "Deletion",
            minor == 0 & (Major > 2 | Major == 0) ~ "Deletion",
            # minor == 0 & Major == 1 ~ "LOH",
            # minor == 0 & Major == 2 ~ "CNLOH",
            # minor == 0 & (Major > 2 | Major == 0) ~ "Deletion",
            minor > 0 & ((Major + minor) >= 3) ~ "Gain",
            TRUE ~ "None"
          )
        )

        # Expand all bins with relevant information
        x_chr = x_chr %>% filter(CNA != "None")

        if(nrow(x_chr) == 0) return(NULL)

        lapply(1:nrow(x_chr),
               function(i) {
                 bins_seq = seq(x_chr$from_bin[i], x_chr$to_bin[i], by = 1) * delta

                 data.frame(
                   chr = chr,
                   from = bins_seq,
                   to = bins_seq + delta,
                   CNA = x_chr$CNA[i]
                 )
               }) %>%
          Reduce(f = bind_rows) %>%
          dplyr::distinct()
      },
      PARAMS = lapply(x$cna$chr %>% unique, list),
      parallel = FALSE,
      filter_errors = TRUE,
      silent = TRUE
    ) %>%
      Reduce(f = bind_rows) %>%
      as_tibble()
  }

  # Check all references
  n_refs = sapply(x, function(x) x$reference_genome) %>% table()

  if(n_refs != length(x)) stop("The input objects have different references")

  # Chop everything
  delta_kb = round(delta/1e3, 0)
  cli::cli_h3("Breaking input segments at {.field {delta_kb} Kb} resolution")

  segments = lapply(x %>% seq_along,
                    function(z)
                      chop_clonal_segments(x[[z]], delta) %>%
                      dplyr::mutate(sample = z)
  ) %>%
    Reduce(f = bind_rows)

  # Count everything
  segments_binned = segments %>%
    dplyr::group_by(chr, from, to, CNA) %>%
    dplyr::summarise(n = n()) %>%
    ungroup()

  # Scale to absolute coordinates
  segments_binned = x[[1]] %>%
    relative_to_absolute_coordinates(segments_binned) %>%
    dplyr::mutate(n = ifelse(CNA == 'Deletion', -1 * n, n))

  # Blank plot
  reference_coordinates = get_reference(x[[1]]$reference_genome)

  low = min(reference_coordinates$from)
  upp = max(reference_coordinates$to)

  segments_plot = ggplot(reference_coordinates) +
    my_ggplot_theme() +
    geom_rect(
      aes(
        xmin = centromerStart,
        xmax = centromerEnd,
        ymin = -Inf,
        ymax = Inf
      ),
      alpha = .3,
      colour = 'gainsboro'
    ) +
    geom_segment(
      data = reference_coordinates,
      aes(
        x = from,
        xend = from,
        y = -Inf,
        yend = Inf
      ),
      size = .1,
      color = 'black',
      linetype = 8
    ) +
    geom_hline(yintercept = 0,
               size = 1,
               colour = 'gainsboro') +
    labs(x = "Chromosome", y = paste0("Cases (out of ", length(x), ')')) +
    scale_x_continuous(
      breaks = c(0, reference_coordinates$from, upp),
      labels = c("", gsub(pattern = 'chr', replacement = '', reference_coordinates$chr), "")
    )

  # Label genome covered
  labels_segments = segments_binned %>%
    dplyr::group_by(CNA) %>%
    dplyr::mutate(Mb = abs(n) * (to-from)) %>%
    dplyr::summarise(Mb = sum(Mb)/1e6) %>%
    dplyr::mutate(label = paste0(CNA, " (", Mb, " Mb)"))

  v_labels_segments = labels_segments$label
  names(v_labels_segments) = labels_segments$CNA

  segments_binned$CNA = v_labels_segments[segments_binned$CNA]

  # Colours
  colours = c(`Deletion` = 'steelblue', `Gain` = 'indianred3')
  names(colours) = v_labels_segments[colours %>% names()]

  # Add info to plot
  segments_plot = segments_plot +
    geom_rect(
      data = segments_binned,
      aes(
        xmin = from,
        xmax = to,
        ymin = 0,
        ymax = n,
        fill = CNA
      ),
      alpha = .8
    ) +
    scale_fill_manual(values = colours) +
    guides(fill = guide_legend('', override.aes = list(alpha = 1)))

  brks_y = segments_binned$n %>% abs %>% max
  brks_step = (brks_y/10) %>% floor
  brks_y = seq(-brks_y, brks_y, ifelse(brks_step < 1, 1, brks_step))
  brks_y_labels = abs(brks_y) %>% sapply(function(r) paste0(r, ' (', round(100 * r/length(x), 0), '%)'))

  segments_plot = segments_plot +
    scale_y_continuous(breaks = brks_y, labels = brks_y_labels)

  # Maybe later we add these...
  # drivers = lapply(x, function(x)
  #   if ("is_driver" %in% colnames(x$mutations))
  #     x$mutations %>% dplyr::filter(is_driver)) %>%
  #   Reduce(f = bind_rows) %>%
  #   dplyr::group_by(chr, from, to, gene, driver_label) %>%
  #   dplyr::summarise(n = n()) %>%
  #   dplyr::ungroup()
  #
  # drivers = relative_to_absolute_coordinates(x[[1]], drivers)
  #
  # segments_plot +
  #   geom_line()

  return(segments_plot)
}

# Circular layout
aux_plot_cohort_CNA_circular = function(x, ...)
{
  Ln = names(L)
  if(is.null(Ln)) {
    Ln = paste0("Sample ", 1:length(L))
    names(L) = Ln

    cli::cli_alert_warning("The input list is un-named, using default naming scheme Sample*")
  }

  KARYO_colors = CNAqc:::get_karyotypes_colors(NULL)

  # Extract calls, and flatten them for plotting
  calls = lapply(Ln,
                 function(s)
                 {
                   W = L[[s]]$cna %>%
                     mutate(
                       label = paste(Major, minor, sep = ':'),
                       CN = minor + Major,
                       sample = s
                     ) %>%
                     select(chr, from, to, label, CN, sample)

                   CNAqc:::relative_to_absolute_coordinates(L[[s]], W)
                 })

  calls_flat =
    suppressWarnings(
      Reduce(
        function(x, y) full_join(x, y, by = c("chr", "from", "to", "label", "CN", "sample")),
        calls) %>%
        mutate(
          label = ifelse(label %in% names(KARYO_colors), label, 'other')
        )
    )

  KARYO_colors = c(KARYO_colors, `other` = 'gray')

  chromosomes = calls_flat$chr %>% unique

  # Reference genome
  reference_genome = CNAqc:::get_reference(L[[1]]$reference_genome) %>% filter(chr %in% chromosomes)
  low = min(reference_genome$from)
  upp = max(reference_genome$to)

  # Default blank genome -- remove labels with label_chr = NA
  bl_genome = suppressMessages(
    CNAqc:::blank_genome(ref = L[[1]]$reference_genome, chromosomes = chromosomes, 
                         label_chr = NA) +
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

