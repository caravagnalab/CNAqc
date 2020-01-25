# data('example_dataset_CNAqc')
# x = CNAqc::init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna, example_dataset_CNAqc$purity)
# CNAqc::plot_segments(x) + ylim(0, 5)
# x = detect_wg_overfragmentation(x,  small_segments_cutoff = 1e6)

#' Title
#'
#' @param x
#' @param binomial_p
#' @param small_segments_cutoff
#' @param alpha
#'
#' @return
#' @export
#'
#' @examples
detect_wg_overfragmentation = function(x,
                                        alpha = 0.05,
                                        binomial_p = 0.2,
                                        small_segments_cutoff = 1e6)
{
  clonal_cna = x$cna %>% filter(CCF == 1)

  n_short = sum(clonal_cna$length < small_segments_cutoff)
  n_long = sum(clonal_cna$length >= small_segments_cutoff)

  cli::cli_alert_info("One-tailed Binomial test, whole genome fragments with cutoff {.value {small_segments_cutoff}} bases: {.value {n_short}} short, {.value {n_long}} long.")

  pvalue = frequentist_test_fragmentation(
    n_short = n_short,
    n_long = n_long,
    chr = 'wg',
    arm = '',
    testable = TRUE,
    p_cutoff_short = binomial_p,
    N_tests = 1,
    alpha = alpha)

  if(pvalue < alpha)
    cli::cli_alert_warning("The whole genome is overfragmented with alpha = {.value {alpha}}; p = {.value {pvalue}}")
  else
    cli::cli_alert_warning("The whole genome is not overfragmented with alpha = {.value {alpha}}; p = {.value {pvalue}}")


  x$wg_fragmentation =
    list(
      pvalue = pvalue,
      alpha = alpha,
      is_overfragmented = pvalue < alpha,
      small_segments_cutoff = small_segments_cutoff)

  return(x)
}


