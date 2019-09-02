#' Analyze the peaks of this
#'
#' @param snvs
#' @param cna
#' @param tumour_purity
#' @param karyotypes
#' @param min_karyotype_size
#' @param adjust
#' @param min_density
#' @param neighlim
#' @param matching_epsilon
#' @param ...
#'
#' @return
#' @export
#'
#' @import tidyverse
#' @import peakPick
#'
#' @examples
analyze_peaks = function(snvs,
                         cna,
                         tumour_purity,
                         karyotypes = c('1:0', '1:1', '2:1', '2:0', '2:2'),
                         min_karyotype_size = 0.05,
                         adjust = 1,
                         min_density = 0.5,
                         neighlim = 1,
                         matching_epsilon = 0.015,
                         ...)
{
  input = prepare_input_data(snvs, cna, tumour_purity)
  snvs = input$snvs
  cna = input$cna

  nsnvs = nrow(snvs)
  ncna = nrow(cna)

  ncnacl = sum(cna$CCF == 1)
  ncnasbcl = sum(cna$CCF < 1)

  pio::pioHdr("CNAqc - CNA Quality Check")
  cat('\n')

  num_mappable = sum(!is.na(snvs$karyotype))
  perc_mappable = round(num_mappable / nsnvs * 100)

  # Karyotypes of interest, and filter for karyotype size
  qc_snvs = snvs %>%
    filter(karyotype %in% karyotypes)

  filtered_qc_snvs = qc_snvs %>%
    group_by(karyotype) %>%
    summarise(n = n(), n_proportion = n() / nsnvs) %>%
    arrange(desc(n)) %>%
    mutate(QC = n_proportion > min_karyotype_size)

  # Re-normalize karyotype size for the one with QC = true
  N_total = sum(filtered_qc_snvs %>% filter(QC) %>% pull(n_proportion))
  filtered_qc_snvs = filtered_qc_snvs %>%
    mutate(norm_prop = ifelse(QC, n_proportion / N_total, NA))

  pio::pioStr(
    "\nMutations.",
    'n =',
    num_mappable,
    "mappable",
    paste0('(~', perc_mappable, '% of input) to'),
    sum(filtered_qc_snvs$QC),
    "karyotypes",
    'with at least',
    round(min_karyotype_size * nsnvs),
    'mutations',
    paste0('(', min_karyotype_size * 100, '% of n)\n\n')
  )

  print(filtered_qc_snvs)

  # Actual data and analysis
  qc_karyotypes = filtered_qc_snvs %>% filter(QC) %>% pull(karyotype)
  qc_snvs = qc_snvs %>% filter(karyotype %in% qc_karyotypes)

  pio::pioStr("\nRunning peak detector with the following parameters\n")
  cat("          Tumour purity =", tumour_purity, '\n')
  cat("KDE bandwith adjustment =", adjust, "and density cutoff =", min_density, '\n')
  cat(" peakPick peak neighlim =", neighlim, "matched with epsilon =", matching_epsilon, '\n')

  detections = NULL
  for (k in karyotypes)
  {
    # For karyotypes that cannot be checked we create a dummy entry
    if (!(k %in% qc_karyotypes))
    {
      detection = list(matching = NULL, plot = ggplot() + geom_blank())
    }
    else
    {
      # The other ones are actually checked
      AB = as.numeric(strsplit(k, ':')[[1]])
      expectation = expected_vaf_peak(AB[1], AB[2], tumour_purity)

      detection = peak_detector(
        snvs = snvs %>% filter(karyotype == k),
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

  plot_summarised = ggplot(assembled_corrections, aes(x = peak, y = offset, color = karyotype)) +
    geom_segment(aes(xend = peak, yend = 0)) +
    geom_point(aes(size = weight), show.legend = FALSE) +
    guides(color = guide_legend('')) +
    theme_minimal() +
    theme(
      legend.position = 'bottom',
      legend.key.height = unit(3, 'mm')) +
    xlim(0, 1) +
    labs(x = 'Peak', y = 'Offset (size as weight)') +
    geom_hline(yintercept = overall_score, color = 'red', linetype = 'dashed') +
    # geom_text(
    #   aes(y = overall_score + 0.05 * overall_score),
    #   x = .9,
    #   color = 'red',
    #   label = round(overall_score, 3)) +
    labs(
      title = "Summary matches",
      subtitle = paste0("Score: ", round(overall_score, 3))
    )

  # Assemble an output image, and add information about each karyotype
  assembled_images =
    ggpubr::ggarrange(
      plotlist =
        append(list(plot_summarised), lapply(detections, function(w) w$plot)),
      nrow = 1,
      ncol = length(detections) + 1
    )


  return(
    list(
      score = overall_score,
      matches = assembled_corrections,
      figure = assembled_images
    )
  )

}
