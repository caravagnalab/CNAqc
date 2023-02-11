#' QC by peak-detection algorithms.
#'
#' @description
#'
#' CNAqc uses peak-detection algorithms to QC data; all leverage the idea that
#' VAFs peaks are known for mutations mapped to a segment with given minor/ major
#' allele copies. CNAqc therefore computes expected peaks, and compares them to
#' peaks detected from data. The theory works with minor modifications for both
#' clonal and subclonal segments. Three distinct algorithms are available, each
#' one working with a different type of copy number segment; all analyses are
#' called by this function, which takes care of running all the suitable algorithms
#' based on the input data.
#'
#' * Simple clonal segments (1:0, 2:0, 1:1, 2:1, 2:2). This QC measures an error for
#' the precision of the current purity estimate, failing a whole sample or a subset
#' of segments the value is over a desired maximum value. The error is determined as
#' a linear combination from the distance between VAF peaks and their theoretical
#' expectation. For this analysis, all mutations mapping across any segment with the
#' same major/minor alleles are pooled. Note that this score can be used to select
#' among alternative copy number solutions, i.e., favouring a solution with lower score.
#' The peaks are determined i) via peak-detection algorithms from the peakPick package,
#' applied to a Gaussian kernel density estimate (gKDE) smooth of the VAF distribution,
#' and ii) via the Bmix Binomial mixture model. Peak-matching (i.e., determining what
#' data peak is closest to the expected peak) has two possible implementations: one
#' matching the closest peaks by euclidean distance, the other ranking peaks from
#' higher to lowr VAFs, and prioritising the former.
#'
#' * Complex clonal segments. The QC procedure for these “general” segments uses only
#' the KDE and, as for simple segments, pools all mutations mapping across any segment
#' with the same major/minor alleles. In this case,t no segment-level or sample-level
#' scores are produced, and complex segment with many matched peaks is likely to be correct.
#'
#' * Subclonal simple segments. The QC procedure for these segments uses the KDE and
#' considers 2 subclones with distinct mixing proportions. Differently from clonal CNAs,
#' however, here the analysis is carried out at the level of each segment, i.e., without
#' pooling segments with the same karyotypes. This makes it possible to use subclonal calls
#' from callers that report segment-specific CCF values, e.g., Battenberg. The model in
#' CNAqc ranks the proposed evolutionary alternatives (linear versus branching) based on
#' the number of matched peaks. A subclonal segment with many matched peaks is likely to be
#' correct.
#'
#' Results from peak-based QC are available via \code{plot_peaks_analysis}, and stored
#' inside the input object.
#'
#' @param x A CNAqc object.
#' @param karyotypes For clonal simple CNAs, the list of segments to test; by default LOH regions (A, AA),
#' diploid regions (AB), and amplification regions (AAB, AABB) are tested, corresponding to
#' \code{'1:0', '1:1', '2:1', '2:0', '2:2'} in "Major:minor" notation.
#' @param min_karyotype_size For clonal simple CNAs, a filter for the segments
#' to test. The segment size is defined based on the number of mutations mapped, this cut is on the proportion
#' relative to the whole set of segments one wishes to analyse (defined by `karyotypes`). For example, by setting
#' `min_karyotype_size = 0.2` one would QC clonal simple CNAs that contain at least 20% of the total mutations.
#' The default of this parameter is `0` (all QCed).
#' @param min_absolute_karyotype_mutations For clonal simple CNAs, as \code{min_karyotype_size} but with a cut
#' measured on absolute mutation counts. For example, by setting `min_absolute_karyotype_mutations = 150` one
#' would QC clonal simple CNAs that contain at least `150` mutations. The default of this parameter is `100`.
#' @param p_binsize_peaks For clonal simple CNAs, peaks detected will be filtered if, in a peak, we map
#' less than \code{p_binsize_peaks * N} mutations. The value \code{N} is obtained couting all mutations that map
#' in all peaks. By default this parameters is `0.005`.
#' @param matching_epsilon Deprecated parameter.
#' @param purity_error For clonal simple CNAs, the purity error tolerance to determine QC pass or fail. This can be
#' set automatically using function \code{auto_tolerance} to optimise the analysis based on a desired rate of false
#' positives matches, as a function of the data coverage and (putative) purity.
#' @param VAF_tolerance For clonal simple CNAs,  a tolerance in comparing bands overlaps which is applied
#' to the raw VAF values.
#' @param n_bootstrap For clonal simple CNAs, the number of times peak detection is bootstrapped (by default 1).
#' This helps sometimes finding peaks that might be visually observable but fail to be detected by the underlying
#' peak-detection heuristics.
#' @param matching_strategy For clonal simple CNAs, if \code{"closest"} the closest peak will be used to
#' match the expected peak. If \code{"rightmost"} peaks are matched prioritizing
#' right to left peaks (the higher-VAF gets matched first); this strategy is more correct
#' in principle but works only if there are no spurious peaks in the estimated
#' density. By default the \code{"closest"} strategy is used.
#' @param kernel_adjust For KDE-based matches the adjust density parameter; see \code{density}. Note that a
#' Gaussian kernel is used by setting (\code{kernel = 'gaussian'}).
#' @param KDE Deprecated parameter.
#' @param starting_state_subclonal_evolution For subclonal simple CNAs, the starting state to determine linear versus
#' branching evolutionary models. By default this is an heterozygous diploid `1:1` state.
#' @param cluster_subclonal_CCF For subclonal segments, should the tool try to merge segments with similar CCF and the same copy number alteration?
#'
#' @return An object of class \code{cnaqc}, modified to hold the results from this analysis. For every type
#' of segment analyzed tables with summary peaks are available in \code{x$peaks_analysis}. The most helpful table
#' is usually the one for simple clonal CNAs `x$peaks_analysis$matches`, which reports several information:
#'
#' - `mutation_multiplicity`, the number of copies of the mutation (i.e., a phasing information);
#' - `peak`, `x`, `y` the expected peak, and the matched peak (`x` and `y`);
#' - `offset`, `weight` and `score`, the factors of the final score;
#' - `QC`, a pass/fail status for the peak.
#'
#' The overall sample-level QC result is available in `x$peaks_analysis$QC`.
#'
#' @seealso \code{auto_tolerance}, \code{plot_peak_analysis} and \code{plot_QC}.
#'
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(mutations = example_dataset_CNAqc$mutations, cna = example_dataset_CNAqc$cna, purity = example_dataset_CNAqc$purity)
#'
#' # Note the run outputs
#' x = analyze_peaks(x)
#'
#' # More precise messages
#' print(x)
#'
#' # The tabulars with summary results per peak and segment
#' print(x$peaks_analysis)
#'
#' # Analysis where simple clonal segments are matched with an alternative algorithm.
#' x = analyze_peaks(x, matching_strategy = "rightmost")
#'
#' print(x)
analyze_peaks = function(x,
                         karyotypes = c('1:0', '1:1', '2:0', '2:1', '2:2'),
                         min_karyotype_size = 0,
                         min_absolute_karyotype_mutations = 100,
                         p_binsize_peaks = 0.005,
                         matching_epsilon = NULL,
                         purity_error = 0.05,
                         VAF_tolerance = 0.015,
                         n_bootstrap = 1,
                         kernel_adjust = 1,
                         matching_strategy = "closest",
                         KDE = TRUE,
                         starting_state_subclonal_evolution = "1:1",
                         cluster_subclonal_CCF = FALSE)
{

  if (!is.null(matching_epsilon)) {
    stop("matching_epsilon is deprecated - using purity_error = ",
         purity_error)
    matching_epsilon = purity_error
  }

  # Check inputs
  stopifnot(inherits(x, "cnaqc"))
  stopifnot(min_karyotype_size >= 0 & min_karyotype_size < 1)
  stopifnot(min_absolute_karyotype_mutations >= 0)
  stopifnot(p_binsize_peaks > 0 & p_binsize_peaks < 1)
  stopifnot(purity_error > 0 & purity_error < 1)
  stopifnot(is.numeric(kernel_adjust))
  stopifnot(matching_strategy %in% c("closest", "rightmost"))

  # Common peaks analysis - they must be in the sample
  cli::cli_h1("Peak analysis: simple CNAs")
  cat("\n")

  x = x %>% analyze_peaks_common(
    karyotypes = karyotypes,
    min_karyotype_size = min_karyotype_size,
    min_absolute_karyotype_mutations = min_absolute_karyotype_mutations,
    p_binsize_peaks = p_binsize_peaks,
    purity_error = purity_error,
    VAF_tolerance = VAF_tolerance,
    n_bootstrap = n_bootstrap,
    kernel_adjust = kernel_adjust
  )

  # Generalised peak analysis
  cli::cli_h1("Peak analysis: complex CNAs")
  cat("\n")

  w = x$n_karyotype[!(x$n_karyotype %>% names() %in% karyotypes)]
  w = w[w > min_absolute_karyotype_mutations]

  if(length(w) > 0)
  {
    cli::cli_alert_info(
        "Karyotypes {.field {names(w)}} with >{.field {min_absolute_karyotype_mutations}} mutation(s). Using epsilon = {.field {purity_error}}."
      )

    x = x %>% analyze_peaks_general(
      n_min = min_absolute_karyotype_mutations,
      epsilon = purity_error,
      kernel_adjust = kernel_adjust,
      n_bootstrap = n_bootstrap
      )

    x$peaks_analysis$general$summary %>%
      print()
  }
  else
    cli::cli_alert_info(
      "No karyotypes with >{.field {min_absolute_karyotype_mutations}} mutation(s). "
    )

  cli::cli_h1("Peak analysis: subclonal CNAs")
  cat("\n")

  # Subclonal CNAs peak analysis
  if(x$n_cna_subclonal > 0)
  {
    x = x %>% analyze_peaks_subclonal(
      n_min = min_absolute_karyotype_mutations,
      epsilon = purity_error,
      kernel_adjust = kernel_adjust,
      n_bootstrap = n_bootstrap,
      starting_state = starting_state_subclonal_evolution,
      cluster_subclonal_CCF = cluster_subclonal_CCF
    )

    if(!is.null(x$peaks_analysis$subclonal))
      x$peaks_analysis$subclonal$summary %>% print()
    else
      cli::cli_alert_info("Subclonal CNAs not analysed with the current parameters.")
  }
  else
    cli::cli_alert_info("No subclonal CNAs in this sample.")

  return(x)
}

# analyze_peaks = function(x,
#                          karyotypes = c('1:0', '1:1', '2:0', '2:1', '2:2'),
#                          min_karyotype_size = 0,
#                          min_absolute_karyotype_mutations = 100,
#                          p_binsize_peaks = 0.005,
#                          matching_epsilon = NULL,
#                          purity_error = 0.05,
#                          VAF_tolerance = 0.015,
#                          n_bootstrap = 1,
#                          kernel_adjust = 1,
#                          matching_strategy = "closest",
#                          KDE = TRUE,
#                          starting_state_subclonal_evolution = "1:1")
# {
#
#   if (!is.null(matching_epsilon)) {
#     stop("matching_epsilon is deprecated - using purity_error = ",
#          purity_error)
#     matching_epsilon = purity_error
#   }
#
#   # Check input
#   stopifnot(inherits(x, "cnaqc"))
#   stopifnot(min_karyotype_size >= 0 & min_karyotype_size < 1)
#   stopifnot(min_absolute_karyotype_mutations >= 0)
#   stopifnot(p_binsize_peaks > 0 & p_binsize_peaks < 1)
#   stopifnot(purity_error > 0 & purity_error < 1)
#   stopifnot(is.numeric(kernel_adjust))
#   stopifnot(matching_strategy %in% c("closest", "rightmost"))
#
#   # Generalised peak analysis
#   cli::cli_h1("Peak analysis: sample-level QC with common karyotypes")
#   cat("\n")
#
#
#   # Karyotypes of interest, and filter for karyotype size
#   qc_mutations = x$mutations %>%
#     dplyr::filter(karyotype %in% karyotypes)
#
#   filtered_qc_mutations = qc_mutations %>%
#     dplyr::group_by(karyotype) %>%
#     dplyr::summarise(
#       n = n(),
#       n_proportion = n() / x$n_mutations,
#       .groups = 'drop'
#     ) %>%
#     dplyr::arrange(desc(n)) %>%
#     dplyr::mutate(QC = (n_proportion > min_karyotype_size) &
#                     n > min_absolute_karyotype_mutations)
#
#   # Re-normalize karyotype size for the ones with QC = true
#   N_total = sum(filtered_qc_mutations %>% filter(QC) %>% pull(n_proportion))
#   filtered_qc_mutations = filtered_qc_mutations %>%
#     dplyr::mutate(norm_prop = ifelse(QC, n_proportion / N_total, NA))
#
#   n_k = sum(filtered_qc_mutations %>% filter(QC) %>% pull(n))
#
#   cli::cli_alert_info(
#     "Karyotypes {.field {karyotypes}}. Matching strategy {.field {matching_strategy}}. KDE = {.field {KDE}}."
#   )
#
#   cli::cli_alert_info(
#     paste0(
#       "Found n = ",
#       n_k,
#       ' mutations in ',
#       paste(
#         filtered_qc_mutations %>% dplyr::filter(QC) %>% dplyr::pull(karyotype),
#         collapse = ', '
#       ),
#       ' (skipping those with n < ',
#       round(min_karyotype_size * x$n_mutations),
#       ' mutations).'
#     )
#   )
#
#   # Special case, nothing to compute
#   if (n_k == 0) {
#     warning("There are non mutations to detect peaks, returning input object.")
#     return(x)
#   }
#
#   # Actual data and analysis
#   qc_karyotypes = filtered_qc_mutations %>% dplyr::filter(QC) %>% dplyr::pull(karyotype)
#   qc_mutations = qc_mutations %>% dplyr::filter(karyotype %in% qc_karyotypes)
#
#   tumour_purity = x$purity
#
#   detections = NULL
#   for (k in karyotypes)
#   {
#     # For karyotypes that cannot be checked we create a dummy entry
#     if (!(k %in% qc_karyotypes))
#     {
#       detection = list(matching = NULL)
#     }
#     else
#     {
#       cli::cli_alert("Fitting {.field {k}} from {.field {x$n_karyotype[k]}} mutations - {.field {n_bootstrap}} bootstrap(s).")
#
#       # The other ones are actually checked
#       AB = as.numeric(strsplit(k, ':')[[1]])
#       expectation = CNAqc:::expected_vaf_peak(AB[1], AB[2], tumour_purity)
#
#       band_matching = delta_vaf_karyo(epsilon_error = purity_error,
#                                       purity = x$purity) %>%
#         dplyr::filter(karyotype == k) %>%
#         pull(delta_vaf)
#
#       # Bootstrap function - w is a fake parameter, the N total
#       # is used meaning that if N==1 no bootstrap is computed
#       pd = function(w, n_boot)
#       {
#         w = qc_mutations %>%
#           dplyr::filter(karyotype == k)
#
#         # Bootstrap only if multiple samples are required
#         data_input = w
#         if (n_boot > 1)
#           data_input = w %>% dplyr::sample_n(size = nrow(w), replace = TRUE)
#
#         run_results = NULL
#
#         if (matching_strategy == "rightmost")
#           run_results = peak_detector(
#             mutations = data_input,
#             expectation = expectation,
#             tumour_purity = tumour_purity,
#             filtered_qc_mutations = filtered_qc_mutations,
#             p = p_binsize_peaks,
#             kernel_adjust = kernel_adjust,
#             matching_epsilon = band_matching,
#             VAF_tolerance = VAF_tolerance,
#             KDE = KDE
#           )
#
#         if (matching_strategy == "closest")
#           run_results = peak_detector_closest_hit_match(
#             mutations = data_input,
#             expectation = expectation,
#             tumour_purity = tumour_purity,
#             filtered_qc_mutations = filtered_qc_mutations,
#             p = p_binsize_peaks,
#             kernel_adjust = kernel_adjust,
#             matching_epsilon = band_matching,
#             VAF_tolerance = VAF_tolerance,
#             KDE = KDE
#           )
#
#         run_results
#       }
#
#       # Sample "n_bootstrap" times via the bootstrap
#       if (n_bootstrap < 0)
#         n_bootstrap = 1
#
#       # For error handling, we switch to easypar
#       # best_peaks = lapply(1:n_bootstrap, pd, n_boot = n_bootstrap)
#       best_peaks = easypar::run(
#         FUN = pd,
#         PARAMS = lapply(1:n_bootstrap, function(w)
#           data.frame(w = w, n_boot = n_bootstrap)),
#         parallel = FALSE,
#         silent = TRUE,
#         filter_errors = TRUE
#       )
#
#       # Runs for this giving ALL errors are intercepted and an error is raised
#       # https://github.com/caravagnalab/CNAqc/issues/10
#       # Without a better example failure a more intelligent handling is difficult to achieve
#       if (length(best_peaks) == 0)
#         stop(
#           "Cannot QC karyotype ",
#           k,
#           " - consider removing that from the input \"karyotypes\" vector."
#         )
#
#       best_score = which.min(sapply(best_peaks, function(w)
#         w$matching$weight %*% w$matching$offset))
#       best_score = best_score[1]
#
#       # Best is returned
#       detection = best_peaks[[best_score]]
#     }
#
#     detections = append(detections, list(detection))
#   }
#   names(detections) = karyotypes
#
#   # =-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   # Compute a linear combination for a GOF of each karyotype
#   # =-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   assembled_corrections = Reduce(bind_rows,
#                                  lapply(detections, function(w)
#                                    w$matching)) %>%
#     mutate(score = weight * offset)
#   assembled_corrections$purity_error = purity_error
#
#   overall_score = sum(assembled_corrections$score)
#
#   # ######   ######   ######   ######   ######   ######   ######
#   # # QC of the analysis, each karyotype independently up to the overall score
#   # ######   ######   ######   ######   ######   ######   ######
#   #
#   # # Function to a border to a plot
#   # qc_plot = function(x, QC)
#   # {
#   #   if (is.na(QC))
#   #     return(x)
#   #   qc = ifelse(QC == "FAIL", "indianred3", 'forestgreen')
#   #
#   #   x +
#   #     theme(title = element_text(color = qc),
#   #           panel.border = element_rect(colour = qc,
#   #                                       fill = NA))
#   # }
#
#   # QC per karyotype ~ kscore is karyotype score (sum, partila of the overall score)
#   score_per_karyotype = assembled_corrections %>%
#     dplyr::group_by(karyotype) %>%
#     dplyr::arrange(desc(counts_per_bin)) %>%
#     dplyr::filter(row_number() == 1) %>%
#     dplyr::select(karyotype, matched) %>%
#     dplyr::mutate(QC = ifelse(matched, "PASS", 'FAIL')) %>%
#     dplyr::select(-matched)
#
#   assembled_corrections = assembled_corrections %>% dplyr::left_join(score_per_karyotype, by = 'karyotype')
#
#   print(assembled_corrections)
#
#   # QC overall score
#   # QC = ifelse(abs(overall_score) < matching_epsilon, "PASS", "FAIL")
#
#   QC = assembled_corrections %>%
#     dplyr::group_by(QC) %>%
#     dplyr::summarise(prop = sum(weight)) %>%
#     dplyr::arrange(desc(prop)) %>%
#     dplyr::filter(dplyr::row_number() == 1) %>%
#     dplyr::pull(QC)
#
#   if (QC == "FAIL")
#     cli::cli_alert_danger(
#       "Peak detection {red('FAIL')} with {.value {red(paste0('r = ', overall_score))}} and tolerance e = {.value {2 * purity_error}}"
#     )
#
#   if (QC == "PASS")
#     cli::cli_alert_success(
#       "Peak detection {green('PASS')} with {.value {green(paste0('r = ', overall_score))}} and tolerance e = {.value {2 * purity_error}}"
#     )
#
#   # Plots assembly
#   # qc_each_plot = pio:::nmfy(score_per_karyotype$karyotype, score_per_karyotype$QC)
#   # na_plots = setdiff(names(detections), score_per_karyotype$karyotype)
#   # qc_each_plot[na_plots] = NA
#   #
#   # plots = lapply(names(detections), function(w)
#   #   qc_plot(detections[[w]]$plot, QC = qc_each_plot[w]))
#   # names(plots) = karyotypes
#
#   x$peaks_analysis = list(
#     score = overall_score,
#     fits = detections,
#     matches = assembled_corrections,
#     matching_strategy = matching_strategy,
#     purity_error = purity_error,
#     # plots = plots,
#     min_karyotype_size = min_karyotype_size,
#     p_binsize_peaks = p_binsize_peaks,
#     # matching_epsilon = matching_epsilon,
#     QC = QC,
#     KDE = KDE
#   )
#
#   # Generalised peak analysis
#   cli::cli_h1("Peak analysis: QC with general karyotypes")
#   cat("\n")
#
#   w = x$n_karyotype[!(x$n_karyotype %>% names() %in% karyotypes)]
#   w = w[w > min_absolute_karyotype_mutations]
#
#   if(length(w) > 0)
#   {
#     cli::cli_alert_info(
#       "Karyotypes {.field {names(w)}} with >{.field {min_absolute_karyotype_mutations}} mutation(s). Using epsilon = {.field {purity_error}}."
#     )
#
#     x = x %>% analyze_peaks_general(
#       n_min = min_absolute_karyotype_mutations,
#       epsilon = purity_error,
#       kernel_adjust = kernel_adjust,
#       n_bootstrap = n_bootstrap
#     )
#
#     x$peaks_analysis$general$summary %>%
#       print()
#   }
#   else
#     cli::cli_alert_info(
#       "No karyotypes with >{.field {min_absolute_karyotype_mutations}} mutation(s). "
#     )
#
#
#   cli::cli_h1("Peak analysis: subclonal CNAs")
#   cat("\n")
#
#   # Subclonal CNAs peak analysis
#   if(x$n_cna_subclonal > 0)
#   {
#     x = x %>% analyze_peaks_subclonal(
#       n_min = min_absolute_karyotype_mutations,
#       epsilon = purity_error,
#       kernel_adjust = kernel_adjust,
#       n_bootstrap = n_bootstrap,
#       starting_state = starting_state_subclonal_evolution
#     )
#
#     if(!is.null(x$peaks_analysis$subclonal))
#       x$peaks_analysis$subclonal$summary %>% print()
#     else
#       cli::cli_alert_info("Subclonal CNAs not analysed with the current parameters.")
#   }
#   else
#     cli::cli_alert_info("No subclonal CNAs in this sample.")
#
#   return(x)
# }
