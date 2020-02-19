
#' Analyze CNA calls by peak detection.
#'
#' @description This function carries out a peak-detection
#' analysis based on KDE and the package \code{peakPick}
#' in order to determine if the mutations that map to a
#' certain karyotype have an allelic frequency that is
#' consistent with the Major and minor alleles of the
#' segment, and the tumour purity. A score is produced
#' as a linear combination of the distance of the actual
#' peak to the expected one, derived with standard equations
#' for CNA analysis.
#'
#' @param x An object of class \code{cnaqc}, created by the \code{init} function.
#' @param karyotypes The list of karyotypes to test. By default LOH regions (A, AA),
#' diploid regions (AB), and amplification regions (AAB, AABB). These correspond to
#' \code{'1:0', '1:1', '2:1', '2:0', '2:2'} in "Major:minor" notation.
#' @param min_karyotype_size Karyotypes are subset by their size (normalized for
#' the number of input mutations). Karyotypes smaller than this value are removed
#' from analysis.
#' @param kernel_adjust KDE adjust density parameter; see \code{density}. A Gaussian
#' kernel is used (\code{kernel = 'gaussian'}).
#' @param p_binsize_peaks Peaks detected will be filtered if, in a peak, we map
#' less than \code{p_binsize_peaks * N} mutations. The value \code{N} is obtained
#' couting all mutations that map in all peaks.
#' @param matching_epsilon Peaks at a location are matched with a range of \code{epsilon * 0.5},
#' left and right. By default (\code{0.025}) a maxium 5\% tolerance is adopted. The status of
#' overall PASS/ FAIL for this analysis is determined with the same tolerance, given the score
#' computed across all karyotypes.
#'
#' @return An object of class \code{cnaqc}, modified to hold the results from this analysis.
#' See the vignette to see how to extract and plot the results.
#' @export
#'
#' @import peakPick
#' @import ggrepel
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna,example_dataset_CNAqc$purity)
#'
#' x = analyze_peaks(x)
#' print(x)
#'
#' print(x$peaks_analysis)
analyze_peaks = function(x,
                         karyotypes = c('1:0', '1:1', '2:1', '2:0', '2:2'),
                         min_karyotype_size = 0.05,
                         p_binsize_peaks = 0.005,
                         matching_epsilon = 0.025,
                         kernel_adjust = 1
                         )
{
  # Check input
  stopifnot(inherits(x, "cnaqc"))
  stopifnot(min_karyotype_size >0 & min_karyotype_size < 1)
  stopifnot(p_binsize_peaks >0 & p_binsize_peaks < 1)
  stopifnot(matching_epsilon >0 & matching_epsilon < 1)
  stopifnot(is.numeric(kernel_adjust))

  # Karyotypes of interest, and filter for karyotype size
  qc_snvs = x$snvs %>%
    dplyr::filter(karyotype %in% karyotypes)

  filtered_qc_snvs = qc_snvs %>%
    dplyr::group_by(karyotype) %>%
    summarise(n = n(), n_proportion = n() / x$n_snvs) %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::mutate(QC = n_proportion > min_karyotype_size)

  # Re-normalize karyotype size for the one with QC = true
  N_total = sum(filtered_qc_snvs %>% filter(QC) %>% pull(n_proportion))
  filtered_qc_snvs = filtered_qc_snvs %>%
    dplyr::mutate(norm_prop = ifelse(QC, n_proportion / N_total, NA))

  n_k = sum(filtered_qc_snvs %>% filter(QC) %>% pull(n))

  cli::cli_alert_info(
    paste0(
      "Requested karyotypes ",
      paste(karyotypes, collapse = ', '),
      "; found n = ",
      n_k,
      ' mutations in ',
      paste(
        filtered_qc_snvs %>% dplyr::filter(QC) %>% dplyr::pull(karyotype),
        collapse = ', '
      ),
      ' (skipping those with n < ',
      round(min_karyotype_size * x$n_snvs),
      ' mutations).'
    )
  )

  # Actual data and analysis
  qc_karyotypes = filtered_qc_snvs %>% filter(QC) %>% pull(karyotype)
  qc_snvs = qc_snvs %>% filter(karyotype %in% qc_karyotypes)

  tumour_purity = x$purity


  detections = NULL
  for (k in karyotypes)
  {
    # For karyotypes that cannot be checked we create a dummy entry
    if (!(k %in% qc_karyotypes))
    {
      detection = list(matching = NULL,
                       plot = CNAqc:::eplot())
    }
    else
    {
      # The other ones are actually checked
      AB = as.numeric(strsplit(k, ':')[[1]])
      expectation = CNAqc:::expected_vaf_peak(AB[1], AB[2], tumour_purity)

      detection = peak_detector(
        snvs = qc_snvs %>% dplyr::filter(karyotype == k),
        expectation = expectation,
        tumour_purity = tumour_purity,
        filtered_qc_snvs = filtered_qc_snvs,
        p = p_binsize_peaks,
        kernel_adjust = kernel_adjust,
        matching_epsilon = matching_epsilon)
    }

    detections = append(detections, list(detection))
  }
  names(detections) = karyotypes

  # =-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # Compute a linear combination for a GOF of each karyotype
  # =-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  assembled_corrections = Reduce(
    bind_rows,
    lapply(detections, function(w) w$matching)
    ) %>%
    mutate(score = weight * offset)

  overall_score = sum(assembled_corrections$score)

  ######   ######   ######   ######   ######   ######   ######
  # QC of the analysis, each karyotype independently up to the overall score
  ######   ######   ######   ######   ######   ######   ######

  # Function to a border to a plot
  qc_plot = function(x, QC)
  {
    if(is.na(QC)) return(x)
    qc = ifelse(QC == "FAIL", "indianred3", 'forestgreen')

    x +
      theme(title = element_text(color = qc),
            panel.border = element_rect(
              colour = qc,
              fill = NA,
              size = 5
            ))
  }

  # QC per karyotype ~ kscore is karyotype score (sum, partila of the overall score)
  score_per_karyotype = assembled_corrections %>%
    dplyr::group_by(karyotype) %>%
    dplyr::summarise(kscore = sum(score)) %>%
    dplyr::mutate(QC = ifelse(abs(kscore) < 2 * matching_epsilon, "PASS", "FAIL"))

  assembled_corrections = assembled_corrections %>% dplyr::left_join(score_per_karyotype, by = 'karyotype')

  print(assembled_corrections)

  # QC overall score
  QC = ifelse(abs(overall_score) < 2 * matching_epsilon, "PASS", "FAIL")

  if(QC == "FAIL")
    cli::cli_alert_danger("Peak detection {red('FAIL')} with {.value {red(paste0('r = ', overall_score))}} and tolerance e = {.value {2 * matching_epsilon}}")
  if(QC == "PASS")
    cli::cli_alert_success("Peak detection {green('PASS')} with {.value {green(paste0('r = ', overall_score))}} and tolerance e = {.value {2 * matching_epsilon}}")

  # Plots assembly
  qc_each_plot = pio:::nmfy(score_per_karyotype$karyotype, score_per_karyotype$QC)
  na_plots = setdiff(names(detections), score_per_karyotype$karyotype)
  qc_each_plot[na_plots] = NA

  plots = lapply(names(detections), function(w) qc_plot(detections[[w]]$plot, QC = qc_each_plot[w]))
  names(plots) = karyotypes

  x$peaks_analysis = list(
    score = overall_score,
    matches = assembled_corrections,
    plots = plots,
    QC = QC
  )

  return(x)
}
