# data('example_dataset_CNAqc')
# x = CNAqc::init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna, example_dataset_CNAqc$purity)
# CNAqc::plot_segments(x) + ylim(0, 5)

#' Title
#'
#' @param x
#' @param alpha
#' @param genome_percentage_cutoff
#' @param minimum_segments_for_testing
#'
#' @return
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc')
#' x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna, example_dataset_CNAqc$purity)
#'
#' x = detect_arm_overfragmentation(x)
#' print(x)
detect_arm_overfragmentation = function(x,
                                    alpha = 0.01,
                                    genome_percentage_cutoff = .2,
                                    minimum_segments_for_testing = 10)
{
  clonal_cna = x$cna %>% filter(CCF == 1)

  # Split calls in arms p and q
  expanded_reference = CNAqc:::expand_reference_chr_to_arms(x)

  # Chromosome length
  L = pio:::nmfy(
    expanded_reference$chr,
    expanded_reference$length)

  # Break segments by arm
  clonal_cna = CNAqc:::split_cna_to_arms(x, CNAqc:::relative_to_absolute_coordinates(x, clonal_cna)) %>%
    mutate(
      L = L[paste0(chr, arm)],
      perc_length = length/L,
      smaller =  perc_length <= genome_percentage_cutoff
    )

  # Extract the counts of segments sizes
  counts = clonal_cna %>%
    group_by(chr, arm, smaller) %>%
    summarise(n_short = n()) %>%
    ungroup()

  # Detect overfragmentation exception - no short segments (will not compute)
  if(all(counts$smaller == FALSE)) {
    cli::cli_alert_warning("No short segments with these parameters.")
    return(x)
  }

  # Detect overfragmentation exception - only short segments
  if(all(counts$smaller == TRUE))
  {
    counts = counts %>%
      spread(smaller, n_short) %>%
      rename(n_short = `TRUE`)

    # Force this in
    counts$n_long = 0
  }
  else
  {
    counts = counts %>%
      spread(smaller, n_short) %>%
      rename(n_short = `TRUE`, n_long = `FALSE`)
  }

  # NAs are 0s
  counts$n_long[is.na(counts$n_long)] = 0
  counts$n_short[is.na(counts$n_short)] = 0

  # Jumos
  counts$jumps = apply(
    counts,
    1,
    function(x) CNAqc:::compute_jumps_segments(clonal_cna, x['chr'], x['arm'])
    )


  # Test only those entries with at least minimum_segments_for_testing segments
  testable = counts %>%
    mutate(total_segments = n_short + n_long) %>%
    mutate(testable = total_segments >= minimum_segments_for_testing) %>%
    arrange(desc(testable), desc(total_segments))

  # Test p-value
  N_tests = sum(testable$testable)

  cli::cli_alert_info(
    "One-tailed Binomial test: {.value {N_tests}} tests, alpha {.value {alpha}}. Short segments: {.value {genome_percentage_cutoff}}% of the reference arm."
  )

  testable = testable %>%
    rowwise() %>%
    mutate(
      p_value = frequentist_test_fragmentation(n_short,
                                         n_long,
                                         chr,
                                         arm,
                                         testable,
                                         p_cutoff_short = genome_percentage_cutoff,
                                         N_tests = N_tests,
                                         alpha = alpha),
      Bonferroni_cutoff = alpha/N_tests,
      significant = ifelse(N_tests == 0, FALSE, p_value < (alpha/N_tests))
    ) %>%
    arrange(p_value) %>%
    ungroup()

  N_sign = sum(testable$significant)

  cli::cli_alert_info(
    "{.value {N_sign}} significantly overfragmented chromosome arms (alpha level {.value {alpha}})."
  )

  x$arm_fragmentation =
    list(
      table = testable,
      alpha = alpha,
      N_tests = N_tests,
      genome_percentage_cutoff = genome_percentage_cutoff,
      minimum_segments_for_testing = minimum_segments_for_testing)

  return(x)
}


