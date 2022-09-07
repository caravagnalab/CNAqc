#' Compute CCF values.
#'
#' @description
#'
#' With this function, CNAqc can compute a CCF per mutation upon “phasing” the multiplicity
#' for every input mutation. Phasing is the task of computing the number of copies of a
#' mutation mapping in a certain copy number segment; this task is a difficult, and can lead to erroneous CCF estimates.
#'
#' CNAqc computes CCFs for simple clonal CNA segments, offering two algorithms to phase
#' mutations directly from VAFs.
#'
#' * Entropy based method. The entropy-based approach will flag mutations for which
#'  we cannot phase multiplicity by VAFs with certainty; the CCFs of these mutations
#'  should be manually controlled and, unless necessary, discarded. To this aid, a QC pass
#'  is assigned with less than a certain percentage of mutations have uncertain CCFs. The model uses the
#'  entropy of a VAF mixture with two Binomial distributions to detect mutations happened
#'  before, and after aneuploidy. Assigning multiplicities is difficult at the crossing
#'  of the two densities where mutations could have multiplicity 1 or 2. If mistaken,
#'  these mutations can determine aritficial peaks in the CCF distribution and compromise
#'  downstream subclonal deconvolution.
#'
#'  * Hard-cut based method. A method is available to compute CCFs regardless of the
#'  entropy. From the 2-class Binomial mixture, CNAqc uses the means of the Binomial
#'  parameters to determine a hard split of the data. Since there are no NA assignments,
#'  the computation is always scored PASS for QC purposes; for this reason this computation is more “rough” than the one based on entropy.
#'
#' Like for other analyses This function creates a field `CCF_estimates` inside
#' the returned object which contains the estimated CCFs.
#'
#' @param x A CNAqc object.
#' @param karyotypes The karyotypes to use, this package supports only clonal simple CNAs.
#' @param cutoff_QC_PASS For the entropy-based method, percentage of mutations that
#' can be not-assigned (\code{NA}) in a karyotype. If the karyotype has more than
#' \code{cutoff_QC_PASS} percentage of non-assigned mutations, then the overall set of CCFs
#' is failed for the karyotype.
#' @param muts_per_karyotype Minimum number of mutations that are required to be mapped to a karyotype
#' in order to compute CCF values (default 25).
#' @param method Either \code{"ENTROPY"} (default) or \code{"ROUGH"}, to reflect the two different algorithms
#' to compute CCF.
#'
#' @seealso Getters function \code{CCF} and \code{plot_CCF}.
#' @return A CNAqc object with CCF values.
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc')
#' x = init(mutations = example_dataset_CNAqc$mutations, cna =example_dataset_CNAqc$cna, purity = example_dataset_CNAqc$purity)
#'
#' x = compute_CCF(x, karyotypes = c('1:0', '1:1', '2:1', '2:0', '2:2'))
#' print(x)
#'
#' # Extract the values with these other functions
#' CCF(x)
#' plot_CCF(x)
compute_CCF = function(x,
                       karyotypes = c('1:0', '1:1', '2:0', '2:1', '2:2'),
                       muts_per_karyotype = 25,
                       cutoff_QC_PASS = 0.1,
                       method = 'ENTROPY')
{
  stopifnot(inherits(x, 'cnaqc'))
  stopifnot(method %in% c('ENTROPY', "ROUGH"))

  if(any(x$n_karyotype <= muts_per_karyotype)) warning("Some karyotypes have fewer than", muts_per_karyotype, 'and will not be analysed.')
  nkaryotypes = x$n_karyotype[x$n_karyotype > muts_per_karyotype]

  karyotypes = intersect(karyotypes, nkaryotypes %>% names)
  stopifnot(
    karyotypes %in% c('1:0', '1:1', '2:1', '2:0', '2:2')
  )

  # Compute mutation multiplicity
  x$CCF_estimates = lapply(
    karyotypes,
    function(k)
    {
      if(k %in% c('1:0', '1:1'))
        return(suppressWarnings(CNAqc:::mutmult_single_copy(x, k)))

      if(method == "ENTROPY")
        return(suppressWarnings(CNAqc:::mutmult_two_copies_entropy(x, k)))
      else
        return(suppressWarnings(CNAqc:::mutmult_two_copies_rough(x, k)))
    })
  names(x$CCF_estimates) = karyotypes

  # Check if there is any null (errors), and remove it
  null_entries = sapply(x$CCF_estimates, function(x) all(is.null(x)))
  x$CCF_estimates = x$CCF_estimates[!null_entries]

  # On extreme cases where there is NO CCF available, we just return x
  if(length(x$CCF_estimates) == 0) {
    x$CCF_estimates = NULL
    return(x)
  }

  # Report some stats
  mutations = lapply(x$CCF_estimates , function(x) x$mutations)
  mutations = Reduce(dplyr::bind_rows, mutations)

  # pioDisp(
  #   mutations %>%
  #     dplyr::group_by(karyotype, mutation_multiplicity) %>%
  #     dplyr::summarise(assignments = n()) %>%
  #     dplyr::ungroup()
  # )

  # QC the findings
  N = mutations %>%
    dplyr::group_by(karyotype) %>%
    dplyr::summarise(N = n()) %>%
    dplyr::ungroup()

  NA_N = mutations %>%
    dplyr::group_by(karyotype, mutation_multiplicity) %>%
    dplyr::summarise(Unknown = n()) %>%
    dplyr::filter(is.na(mutation_multiplicity)) %>%
    dplyr::select(-mutation_multiplicity) %>%
    dplyr::ungroup()

  QC_table = N %>%
    dplyr::full_join(NA_N, by = 'karyotype') %>%
    dplyr::mutate(
      Unknown = ifelse(is.na(Unknown), 0, Unknown),
      p_Unkown = Unknown/N,
      QC = ifelse(p_Unkown < cutoff_QC_PASS, "PASS", "FAIL"),
      method = method
    )

  if(any(QC_table$QC == "FAIL"))
  {
    cat('\n')
    cli::cli_h2("Summary CCF assignments. (>{.field {cutoff_QC_PASS*100}%} NAs: not assignable with confidence)")
    print(QC_table)
  }

  for(k in QC_table$karyotype)
    x$CCF_estimates[[k]]$QC_table = QC_table %>% dplyr::filter(karyotype == !!k)

  x
}



