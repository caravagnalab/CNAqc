
analyze_peaks = function(x,
                         karyotypes = c('1:0', '1:1', '2:1', '2:0', '2:2'),
                         min_karyotype_size = 0.05,
                         adjust = 1,
                         min_density = 0.5,
                         neighlim = 1,
                         matching_epsilon = 0.015,
                         ...)
{
  stopifnot(inherits(x, "cnaqc"))

  pio::pioHdr("QC analysis with peaks detection")

  cat('\n')
  print(x)

  # Karyotypes of interest, and filter for karyotype size
  qc_snvs = x$snvs %>%
    filter(karyotype %in% karyotypes)

  filtered_qc_snvs = qc_snvs %>%
    group_by(karyotype) %>%
    summarise(n = n(), n_proportion = n() / x$n_snvs) %>%
    arrange(desc(n)) %>%
    mutate(QC = n_proportion > min_karyotype_size)

  # Re-normalize karyotype size for the one with QC = true
  N_total = sum(filtered_qc_snvs %>% filter(QC) %>% pull(n_proportion))
  filtered_qc_snvs = filtered_qc_snvs %>%
    mutate(norm_prop = ifelse(QC, n_proportion / N_total, NA))

  n_k = sum(filtered_qc_snvs %>% filter(QC) %>% pull(n))


  pio::pioTit(
    "Analysing",
    "~ karyotypes", paste(karyotypes, collapse = ', '),
    "~ ", n_k, 'mutations',
    "~ min. k =", round(min_karyotype_size * x$n_snvs))

  print(filtered_qc_snvs)

  # Actual data and analysis
  qc_karyotypes = filtered_qc_snvs %>% filter(QC) %>% pull(karyotype)
  qc_snvs = qc_snvs %>% filter(karyotype %in% qc_karyotypes)

  tumour_purity = x$purity

  pio::pioStr(
    "\nPeak detector",
    paste0(
      'p = ', tumour_purity,
      ' ~ KDE a = ', adjust, ' c = ', min_density,
      ' ~ peakPick n = ', neighlim, ' epsilon = ', matching_epsilon
    ), suffix = '\n')

  e = function(){
    ggplot() +
      geom_blank() +
      theme(
        plot.background = element_rect(fill = "lightgray"),
        panel.background = element_rect(fill = "lightgray"),
        )
  }

  detections = NULL
  for (k in karyotypes)
  {
    # For karyotypes that cannot be checked we create a dummy entry
    if (!(k %in% qc_karyotypes))
    {
      detection = list(matching = NULL,
                       plot = e())
    }
    else
    {
      # The other ones are actually checked
      AB = as.numeric(strsplit(k, ':')[[1]])
      expectation = expected_vaf_peak(AB[1], AB[2], tumour_purity)

      detection = peak_detector(
        snvs = qc_snvs %>% filter(karyotype == k),
        expectation,
        tumour_purity,
        filtered_qc_snvs,
        adjust = adjust,
        min_density = min_density,
        neighlim = neighlim,
        matching_epsilon = matching_epsilon,
        ...
      )

    }

    detections = append(detections, list(detection))
  }
  names(detections) = karyotypes

  # =-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # Compute a linear combination for a GOF of each karyotype
  # =-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  assembled_corrections = lapply(detections, function(w)
    w$matching)
  assembled_corrections = Reduce(bind_rows, assembled_corrections) %>%
    mutate(score = weight * offset)

  overall_score = sum(assembled_corrections$score)

  pio::pioTit("Results table")
  print(assembled_corrections)

  pio::pioStr("\nFit score:", overall_score, '\n')


  #
  #   offset = 0
  #   for(k in qc_karyotypes)
  #   {
  #     sum_offsets = sum(abs(detections[[k]]$matching$offset))
  #     weight = filtered_qc_snvs %>% filter(karyotype == k) %>% pull(norm_prop)
  #
  #     filtered_qc_snvs = filtered_qc_snvs %>%
  #       mutate(
  #         adjustment = ifelse(karyotype == k, sum_offsets * weight, adjustment)
  #       )
  #
  #     offset  = offset + sum_offsets * weight
  #   }

  # plot_summarised = ggplot(assembled_corrections, aes(x = peak, y = offset, color = karyotype)) +
  #   geom_segment(aes(xend = peak, yend = 0)) +
  #   geom_point(aes(size = weight), show.legend = FALSE) +
  #   guides(color = guide_legend('')) +
  #   theme_minimal() +
  #   theme(
  #     legend.position = 'bottom',
  #     legend.key.height = unit(3, 'mm')) +
  #   xlim(0, 1) +
  #   labs(x = 'Peak', y = 'Offset (size as weight)') +
  #   geom_hline(yintercept = overall_score, color = 'red', linetype = 'dashed') +
  #   # geom_text(
  #   #   aes(y = overall_score + 0.05 * overall_score),
  #   #   x = .9,
  #   #   color = 'red',
  #   #   label = round(overall_score, 3)) +
  #   labs(
  #     title = "Summary matches",
  #     subtitle = paste0("Score: ", round(overall_score, 3))
  #   )

  plots = lapply(detections, function(w) w$plot)
  names(plots) = karyotypes

  # Assemble an output image, and add information about each karyotype
  # assembled_images =
  #   ggpubr::ggarrange(
  #     plotlist =
  #       append(list(plot_summarised), lapply(detections, function(w) w$plot)),
  #     nrow = 1,
  #     ncol = length(detections) + 1
  #   )


  x$peaks_analysis = list(
    score = overall_score,
    matches = assembled_corrections,
    plots = plots
  )

  return(x)

}
