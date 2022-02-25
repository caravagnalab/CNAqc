#' Randomly subsample a dataset's mutations.
#'
#' @description
#'
#' This functions subsample mutations with uniform probability, retaining
#' all the copy number calls. If data contains driver mutation annotations,
#' these can be forced to remain.
#'
#' @param x An object of class \code{cnaqc}, created by the \code{init} function.
#' @param N The maximum number of mutations to retain.
#' @param keep_drivers If \code{TRUE}, it retains drivers annotated in the data.
#'
#' @param x An object of class \code{cnaqc}, created by the \code{init} function.
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna,example_dataset_CNAqc$purity)
#'
#' subsample(x, N = 100)
#' subsample(x, N = 1000)
subsample = function(x, N = 15000, keep_drivers = TRUE)
{
  if (x$n_snvs < N)
    return(x)

  xx = x$snvs %>% dplyr::sample_n(N)

  xx_d = data.frame(stringsAsFactors = FALSE)
  if (keep_drivers & 'is_driver' %in% colnames(x$snvs))
    xx_d = x$snvs %>% dplyr::filter(is_driver)

  CNAqc::init(
    snvs = dplyr::bind_rows(xx_d, xx),
    cna = x$cna %>% dplyr::select(-segment_id,-n,-CCF),
    purity = x$purity
  )
}

#' Subset calls by segments' karyotype.
#'
#' @description
#'
#' For an object already created, it subsets the calls to those that involve
#' segments with a certain karyotype.
#'
#' @param x An object of class \code{cnaqc}, created by the \code{init} function.
#' @param karyotype A list of karyotype ids in \code{"Major:minor"} notation
#' (e.g., \code{"1:1", "2,1", ...}) that will be retained.
#'
#' @return An object of class \code{cnaqc}, created by the \code{init} function.
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna,example_dataset_CNAqc$purity)
#'
#' subset_by_segment_karyotype(x, '2:2')
#' subset_by_segment_karyotype(x, '2:1')
subset_by_segment_karyotype = function(x, karyotypes)
{
  if (!inherits(x, "cnaqc"))
    stop("Not a CNAqc object in input.")

  cna_calls = x$cna %>%
    dplyr::mutate(k = paste0(Major, ':', minor)) %>%
    dplyr::filter(k %in% karyotypes)

  if (nrow(cna_calls) == 0)
    stop("There are no calls with these karyotypes, cannot subset.")

  return(
    CNAqc::init(
      snvs = x$snvs,
      cna = cna_calls %>% dplyr::select(-segment_id,-n,-CCF),
      purity = x$purity,
      ref = x$reference_genome
    )
  )
}


#' Subset calls by segments total copy number
#'
#' @description
#'
#' For an object already created, it subsets the calls to those that involve
#' segments with a certain total copy number (defined as Major plus minor).
#'
#' @param x An object of class \code{cnaqc}, created by the \code{init} function.
#' @param total A list of copy number integers to find segments to retain only if they
#' have a certain total copy state.
#'
#' @return An object of class \code{cnaqc}, created by the \code{init} function.
#'
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna,example_dataset_CNAqc$purity)
#'
subset_by_segment_totalcn = function(x, totalcn)
{
  if (!inherits(x, "cnaqc"))
    stop("Not a CNAqc object in input.")

  cna_calls = x$cna %>%
    dplyr::mutate(total_cn = Major + minor) %>%
    dplyr::filter(total_cn %in% totalcn)

  if (nrow(cna_calls) == 0)
    stop("There are no calls with these total copy states, cannot subset.")

  return(
    CNAqc::init(
      snvs = x$snvs,
      cna = cna_calls,
      purity = x$purity,
      ref = x$reference_genome
    )
  )
}

subset_by_segment_minmutations = function(x, totalcn)
{
  if (!inherits(x, "cnaqc"))
    stop("Not a CNAqc object in input.")

  cna_calls = x$cna %>%
    dplyr::mutate(total_cn = Major + minor) %>%
    dplyr::filter(total_cn %in% totalcn)

  if (nrow(cna_calls) == 0)
    stop("There are no calls with these total copy states, cannot subset.")

  return(
    CNAqc::init(
      snvs = x$snvs,
      cna = cna_calls,
      purity = x$purity,
      ref = x$reference_genome
    )
  )

}

#' Title
#'
#' @param x
#' @param min_target_CCF
#'
#' @return
#' @export
#'
#' @examples
#'
#' x = init(snvs = example_dataset_CNAqc$snvs, cna = example_dataset_CNAqc$cna, purity = example_dataset_CNAqc$purity)
#' print(x)
#'
#' plot_data_histogram(x)
#'
#' plot_data_histogram(subset_by_minimum_CCF(x))
subset_by_minimum_CCF = function(x, min_target_CCF = 0.1)
{
  used_karyotypes = x$snvs$karyotype %>% unique

  minor_copies = sapply(strsplit(used_karyotypes, split = ':'), function(x)
    x[2]) %>% as.numeric()
  major_copies = sapply(strsplit(used_karyotypes, split = ':'), function(x)
    x[1]) %>% as.numeric()

  purity = x$purity
  mut.allele = 1

  cutoffs = sapply(seq_along(used_karyotypes), function(i) {
    vaf_from_ccf(
      ccf = min_target_CCF,
      m = minor_copies[i],
      M = major_copies[i],
      p = purity,
      mut.allele = mut.allele
    )
  })
  names(cutoffs) = used_karyotypes

  cuts_table = data.frame(
    karyotype = used_karyotypes,
    min_target_CCF = min_target_CCF,
    VAF_cutoff = cutoffs,
    stringsAsFactors = FALSE
  )
  cuts_table$n = x$n_karyotype[cuts_table$karyotype]

  observed_minima = x$snv %>%
    dplyr::group_by(karyotype) %>%
    dplyr::summarise(VAF_minimum = min(VAF), .groups = 'drop') %>%
    dplyr::full_join(cuts_table, by = 'karyotype') %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::select(karyotype, n, VAF_cutoff, VAF_minimum) %>%
    dplyr::mutate(any_to_filter = ifelse(VAF_minimum < VAF_cutoff, "Yes", "No"))

  cli::cli_h1("Cutoffs table for minimum CCF {.field {min_target_CCF}}")
  print(observed_minima)

  # Reports
  will_not_affect = observed_minima %>% dplyr::filter(VAF_cutoff < VAF_minimum)

  if (nrow(will_not_affect) > 0) {
    cat("\n")
    cli::cli_alert_warning(
      "Some karyotypes will not be affected by the desired cut: {.field {will_not_affect$karyotype}}"
    )
  }

  cat("\n")
  subset_data = x$snvs %>%
    dplyr::mutate(remove = VAF <  cutoffs[karyotype]) %>%
    dplyr::filter(!remove)

  init(
    snvs = subset_data,
    cna = x$cna %>% dplyr::select(-segment_id,-n,-CCF),
    purity = x$purity,
    ref = x$reference_genome
  )
}

#' Retain only SNVs.
#'
#' @description It removes all non-SNVs mutations, re-creaging a new CNAqc
#' dataset (all analyses are lost).
#'
#' @param x
#' @param ref_nucleotides What reference alleles to use, default \code{c("A", "C", "T", "G")}.
#' @param alt_nucleotides What alternative alleles to use, default \code{c("A", "C", "T", "G")}.
#'
#' @return A new CNAqc dataset created with \code{init}.
#' @export
#'
#' @examples
#' x = init(snvs = example_dataset_CNAqc$snvs, cna = example_dataset_CNAqc$cna, purity = example_dataset_CNAqc$purity)
#' print(x)
#'
#' # no change - they are already SNVs
#' print(subset_snvs(x))
subset_snvs = function(x,
                       ref_nucleotides = c("A", "C", "T", "G"),
                       alt_nucleotides = c("A", "C", "T", "G"))
{
  subset_snvs = x$snvs %>%
    dplyr::filter(ref %in% ref_nucleotides, alt %in% alt_nucleotides)

  removed_snvs = x$snvs %>%
    dplyr::filter(!(ref %in% ref_nucleotides) |
                    !(alt %in% alt_nucleotides))

  if (nrow(removed_snvs) > 0)
  {
    cli::cli_alert_success(
      "Retained {.field {nrow(subset_snvs)}} SNVs, removed {.field {nrow(removed_snvs)}}."
    )

    x = init(
      snvs = subset_snvs,
      cna = x$cna %>% dplyr::select(-segment_id,-n,-CCF),
      purity = x$purity,
      ref = x$reference_genome
    )
  }

  return(x)
}

# load('/Volumes/Data/Dropbox/160916_HMFregCPCT_FR10302749_FR12251983_CPCT02020357_cnaqc_object.RData')
# plot_data_histogram(x)
# w = init(x$snvs, x$cna, x$purity, x$reference_genome)
# wsb = subset_by_minimum_CCF(w, 0.1)
# wsbsnv = subset_snvs(wsb)
#
# ggpubr::ggarrange(
#   plot_data_histogram(x),
#   plot_data_histogram(wsb),
#   plot_data_histogram(wsbsnv),
#   ncol = 3
#
#
# )


#' Split a dataset by chromosome.
#'
#' @description Split a CNAqc object by chromosome, returning a named list
#' with data splits by chromosomes.
#'
#' @param x An object created by CNAqc.
#' @param chromosomes A list of chromosomes to retain.
#'
#' @return A named list of CNAqc objects.
#' @export
#'
#' @examples
#'
#' data("example_PCAWG", package = 'CNAqc')
#'
#' split_by_chromosome(example_PCAWG)
split_by_chromosome = function(x,
                               chromosomes = paste0('chr', c(1:22, 'X', 'Y')),)
{
  stopifnot(inherits(x, 'cnaqc'))

  objs = NULL
  nm = NULL

  for(chr in chromosomes)
  {
    clonal_mutations = x$snvs %>% filter(chr == !!chr)

    if((clonal_mutations %>% nrow()) == 0) next

    cli::cli_h3(chr)
    cat("\n")

    subclonal_mutations = x$cna_subclonal %>% filter(chr == !!chr) %>% pull(mutations) %>% Reduce(f = bind_rows)
    cnas = x$cna %>% filter(chr == !!chr)

    cnaqc_obj = init(
      clonal_mutations %>% bind_rows(subclonal_mutations),
      cna = cnas,
      purity = x$purity,
      ref = x$reference_genome
    )

    objs = append(objs, list(cnaqc_obj))
    nm = c(nm, chr)
  }

  names(objs) = nm

  return(objs)
}
