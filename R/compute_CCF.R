#' Title
#'
#' @param x
#' @param karyotypes
#' @param entropy_quantile
#'
#' @return
#' @export
#'
#' @examples
compute_CCF = function(x, karyotypes = c('2:1', '2:0', '2:2'), entropy_quantile = .9)
{
  stopifnot(
    karyotypes %in% c('2:1', '2:0', '2:2')
  )

  # Compute mutation multiplicity
  x$CCF_estimates = lapply(karyotypes, mutation_multiplicity_entropy, x = x, entropy_quantile = entropy_quantile)
  names(x$CCF_estimates) = karyotypes

  # Report some stats
  mutations = lapply(x$CCF_estimates , function(x) x$mutations)
  mutations = Reduce(bind_rows, mutations)

  pioTit("Summary CCF assignments")
  pioDisp(
    mutations %>%
      group_by(karyotype, mutation_multiplicity) %>%
      summarise(assignments = n()) %>%
      ungroup()
  )
  cat("Note: NA ~ mutations not confidently assignable at q =", entropy_quantile, '\n')

  x
}

# data('example_dataset_CNAqc')
# x = CNAqc::init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna, example_dataset_CNAqc$purity)
# cc = compute_CCF(x, karyotypes = '2:1')



