

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
#' @param min_absolute_karyotype_mutations As \code{min_karyotype_size} this imposes a
#' minimum number of mutations on a karyotype to be analysed (in absolute number of
#' mutations). Karyotypes smaller than this value are removed from analysis.
#' @param kernel_adjust KDE adjust density parameter; see \code{density}. A Gaussian
#' kernel is used (\code{kernel = 'gaussian'}).
#' @param p_binsize_peaks Peaks detected will be filtered if, in a peak, we map
#' less than \code{p_binsize_peaks * N} mutations. The value \code{N} is obtained
#' couting all mutations that map in all peaks.
#' @param matching_epsilon Peaks at a location are matched with a range of \code{epsilon * 0.5},
#' left and right. By default (\code{0.025}) a maxium 5\% tolerance is adopted. The status of
#' overall PASS/ FAIL for this analysis is determined with the same tolerance, given the score
#' computed across all karyotypes.
#' @param n_bootstrap Number of times peak detection is bootstrapped (this helps sometimes
#' finding peaks that might be visually observable but fail to be detected by the underlying
#' peak detection heuristic). The default of this parameter is 1, meaning that no bootstrap
#' is performed.
#' @param matching_strategy If \code{"closest"}, the closest peak will be used to
#' match the expected peak. If \code{"rightmost"}, peaks are matched prioritizing
#' right to left peaks (the higher get matched first); this strategy is more correct
#' in principle but works only if there are no spurious peaks in the estimated
#' density. By defaule the first strategy is used.
#'
#' @return An object of class \code{cnaqc}, modified to hold the results from this analysis.
#' See the vignette to see how to extract and plot the results.
#' @export
#'
#' @import peakPick
#' @import ggrepel
#' @import BMix
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna,example_dataset_CNAqc$purity)
#'
#' x = analyze_peaks(x)
#' print(x)
#'
#' print(x$peaks_analysis)
#'
#'
#' x = analyze_peaks(x, matching_strategy = "rightmost")
#' print(x)
analyze_peaks = function(x,
                         karyotypes = c('1:0', '1:1', '2:1', '2:0', '2:2'),
                         min_karyotype_size = 0,
                         min_absolute_karyotype_mutations = 10,
                         p_binsize_peaks = 0.005,
                         matching_epsilon = 0.03,
                         n_bootstrap = 1,
                         kernel_adjust = 1,
                         matching_strategy = "closest")
{
  # Check input
  stopifnot(inherits(x, "cnaqc"))
  stopifnot(min_karyotype_size >= 0 & min_karyotype_size < 1)
  stopifnot(min_absolute_karyotype_mutations >= 0)
  stopifnot(p_binsize_peaks > 0 & p_binsize_peaks < 1)
  stopifnot(matching_epsilon > 0 & matching_epsilon < 1)
  stopifnot(is.numeric(kernel_adjust))
  stopifnot(matching_strategy %in% c("closest", "rightmost"))

  # Karyotypes of interest, and filter for karyotype size
  qc_snvs = x$snvs %>%
    dplyr::filter(karyotype %in% karyotypes)

  filtered_qc_snvs = qc_snvs %>%
    dplyr::group_by(karyotype) %>%
    dplyr::summarise(n = n(), n_proportion = n() / x$n_snvs, .groups = 'drop') %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::mutate(QC = (n_proportion > min_karyotype_size) & n > min_absolute_karyotype_mutations)

  # Re-normalize karyotype size for the ones with QC = true
  N_total = sum(filtered_qc_snvs %>% filter(QC) %>% pull(n_proportion))
  filtered_qc_snvs = filtered_qc_snvs %>%
    dplyr::mutate(norm_prop = ifelse(QC, n_proportion / N_total, NA))

  n_k = sum(filtered_qc_snvs %>% filter(QC) %>% pull(n))

  # cli::cli_alert_info(
  #   paste0(
  #     "Requested karyotypes ",
  #     paste(karyotypes, collapse = ', '),
  #     "; found n = ",
  #     n_k,
  #     ' mutations in ',
  #     paste(
  #       filtered_qc_snvs %>% dplyr::filter(QC) %>% dplyr::pull(karyotype),
  #       collapse = ', '
  #     ),
  #     ' (skipping those with n < ',
  #     round(min_karyotype_size * x$n_snvs),
  #     ' mutations).'
  #   )
  # )

  cli::cli_alert_info("Requested karyotypes {.field {karyotypes}}. Matching strategy {.field {matching_strategy}}.")

  cli::cli_alert_info(
    paste0(
      "Found n = ",
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

  # Special case, nothing to compute
  if(n_k == 0) {
    warning("There are non mutations to detect peaks, returning input object.")
    return(x)
  }

  # Actual data and analysis
  qc_karyotypes = filtered_qc_snvs %>% dplyr::filter(QC) %>% dplyr::pull(karyotype)
  qc_snvs = qc_snvs %>% dplyr::filter(karyotype %in% qc_karyotypes)

  tumour_purity = x$purity

  detections = NULL
  for (k in karyotypes)
  {
    # For karyotypes that cannot be checked we create a dummy entry
    if (!(k %in% qc_karyotypes))
    {
      detection = list(matching = NULL)
    }
    else
    {
      # The other ones are actually checked
      AB = as.numeric(strsplit(k, ':')[[1]])
      expectation = CNAqc:::expected_vaf_peak(AB[1], AB[2], tumour_purity)

      # Bootstrap function - w is a fake parameter, the N total
      # is used meaning that if N==1 no bootstrap is computed
      pd = function(w, n_boot)
      {
        w = qc_snvs %>%
          dplyr::filter(karyotype == k)

        # Bootstrap only if multiple samples are required
        data_input = w
        if(n_boot > 1)
          data_input = w %>% dplyr::sample_n(size = nrow(w), replace = TRUE)

        run_results = NULL

        if(matching_strategy == "rightmost")
          run_results = CNAqc:::peak_detector(
            snvs = data_input,
            expectation = expectation,
            tumour_purity = tumour_purity,
            filtered_qc_snvs = filtered_qc_snvs,
            p = p_binsize_peaks,
            kernel_adjust = kernel_adjust,
            matching_epsilon = matching_epsilon
          )

        if(matching_strategy == "closest")
          run_results = CNAqc:::peak_detector_closest_hit_match(
            snvs = data_input,
            expectation = expectation,
            tumour_purity = tumour_purity,
            filtered_qc_snvs = filtered_qc_snvs,
            p = p_binsize_peaks,
            kernel_adjust = kernel_adjust,
            matching_epsilon = matching_epsilon
          )

        run_results
      }

      # Sample "n_bootstrap" times via the bootstrap
      if(n_bootstrap < 0) n_bootstrap = 1

      best_peaks = lapply(1:n_bootstrap, pd, n_boot = n_bootstrap)
      best_score = which.min(sapply(best_peaks, function(w)
        w$matching$weight %*% w$matching$offset))
      best_score = best_score[1]

      # Best is returned
      detection = best_peaks[[best_score]]
    }

    detections = append(detections, list(detection))
  }
  names(detections) = karyotypes

  # =-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # Compute a linear combination for a GOF of each karyotype
  # =-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  assembled_corrections = Reduce(bind_rows,
                                 lapply(detections, function(w)
                                   w$matching)) %>%
    mutate(score = weight * offset)

  overall_score = sum(assembled_corrections$score)

  # ######   ######   ######   ######   ######   ######   ######
  # # QC of the analysis, each karyotype independently up to the overall score
  # ######   ######   ######   ######   ######   ######   ######
  #
  # # Function to a border to a plot
  # qc_plot = function(x, QC)
  # {
  #   if (is.na(QC))
  #     return(x)
  #   qc = ifelse(QC == "FAIL", "indianred3", 'forestgreen')
  #
  #   x +
  #     theme(title = element_text(color = qc),
  #           panel.border = element_rect(colour = qc,
  #                                       fill = NA))
  # }

  # QC per karyotype ~ kscore is karyotype score (sum, partila of the overall score)
  score_per_karyotype = assembled_corrections %>%
    dplyr::group_by(karyotype) %>%
    dplyr::arrange(desc(counts_per_bin)) %>%
    dplyr::filter(row_number() == 1) %>%
    dplyr::select(karyotype, matched) %>%
    dplyr::mutate(QC = ifelse(matched, "PASS", 'FAIL')) %>%
    dplyr::select(-matched)

  assembled_corrections = assembled_corrections %>% dplyr::left_join(score_per_karyotype, by = 'karyotype')

  print(assembled_corrections)

  # QC overall score
  # QC = ifelse(abs(overall_score) < matching_epsilon, "PASS", "FAIL")
  
  QC = assembled_corrections %>% 
    dplyr::group_by(QC) %>% 
    dplyr::summarise(prop = sum(weight)) %>% 
    dplyr::arrange(desc(prop)) %>% 
    dplyr::filter(dplyr::row_number() == 1) %>% 
    dplyr::pull(QC)
  
  if (QC == "FAIL")
    cli::cli_alert_danger(
      "Peak detection {red('FAIL')} with {.value {red(paste0('r = ', overall_score))}} and tolerance e = {.value {2 * matching_epsilon}}"
    )
  
  if (QC == "PASS")
    cli::cli_alert_success(
      "Peak detection {green('PASS')} with {.value {green(paste0('r = ', overall_score))}} and tolerance e = {.value {2 * matching_epsilon}}"
    )

  # Plots assembly
  # qc_each_plot = pio:::nmfy(score_per_karyotype$karyotype, score_per_karyotype$QC)
  # na_plots = setdiff(names(detections), score_per_karyotype$karyotype)
  # qc_each_plot[na_plots] = NA
  #
  # plots = lapply(names(detections), function(w)
  #   qc_plot(detections[[w]]$plot, QC = qc_each_plot[w]))
  # names(plots) = karyotypes

  x$peaks_analysis = list(
    score = overall_score,
    fits = detections,
    matches = assembled_corrections,
    matching_strategy = matching_strategy,
    # plots = plots,
    min_karyotype_size = min_karyotype_size,
    p_binsize_peaks = p_binsize_peaks,
    matching_epsilon = matching_epsilon,
    QC = QC
  )

  return(x)
}
