#' Randomly subsample mutations.
#'
#' @description
#'
#' This functions randoly subsample mutations, retaining
#' all the simple clonal CNAs; subclonal CNAs are dropped.
#' If data contains driver mutation annotations, these can be forced to remain.
#'
#' @param x A CNAqc object.
#' @param N The maximum number of mutations to retain.
#' @param keep_drivers If \code{TRUE}, it retains drivers annotated in the data.
#'
#' @param x A new CNAqc object with subset data.
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(mutations = example_dataset_CNAqc$mutations, cna = example_dataset_CNAqc$cna, purity = example_dataset_CNAqc$purity)
#'
#' # Example runs
#' subsample(x, N = 100)
#' subsample(x, N = 1000)
subsample = function(x, N = 15000, keep_drivers = TRUE)
{
  if (x$n_mutations < N)
    return(x)

  xx = x$mutations %>% dplyr::sample_n(N)

  xx_d = data.frame(stringsAsFactors = FALSE)
  if (keep_drivers & 'is_driver' %in% colnames(x$mutations))
    xx_d = x$mutations %>% dplyr::filter(is_driver)

  init(
    mutations = dplyr::bind_rows(xx_d, xx),
    cna = x$cna %>% dplyr::select(-segment_id,-n,-CCF),
    purity = x$purity,
    sample = x$sample
  )
}

#' Subset by clonal segments.
#'
#' @description
#'
#' Retains only a subset of clonal segments, e.g., all `2:1` segments.
#'
#' @param x A CNAqc object.
#' @param karyotypes A list of karyotype ids in \code{"Major:minor"} notation
#' (e.g., \code{"1:1", "2,1", ...}) that will be retained.
#'
#' @param x A new CNAqc object with subset data.
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(mutations = example_dataset_CNAqc$mutations, cna = example_dataset_CNAqc$cna, purity = example_dataset_CNAqc$purity)
#'
#' subset_by_segment_karyotype(x, '2:2')
#' subset_by_segment_karyotype(x, '2:1')
subset_by_segment_karyotype = function(x, karyotypes)
{
  if (!inherits(x, "cnaqc"))
    stop("Not a CNAqc object in input.")

  cna_calls = x %>% CNA() %>%
    dplyr::mutate(k = paste0(Major, ':', minor)) %>%
    dplyr::filter(k %in% karyotypes)

  if (nrow(cna_calls) == 0)
    stop("There are no calls with these karyotypes, cannot subset.")

  to_remove = c("segment_id", "n", "CCF") %in% colnames(cna_calls)
  to_remove = c("segment_id", "n", "CCF")[to_remove]
  cna_calls %>% dplyr::select(-to_remove)

  return(
    init(
      mutations = x %>% Mutations(),
      cna = cna_calls %>% dplyr::select(to_remove),
      purity = x$purity,
      ref = x$reference_genome,
      sample = x$sample
    )
  )
}


#' Subset clonal simple segments total copy number.
#'
#' @description
#'
#' Retains only a subset of clonal simple segments, e.g., all segments with
#' total ploidy 2.
#'
#' @param x A CNAqc object.
#' @param totalcn The total (Major + minor) copy numbers to filter.
#'
#' @return x A new CNAqc object with subset data.
#'
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(mutations = example_dataset_CNAqc$mutations, cna = example_dataset_CNAqc$cna, purity = example_dataset_CNAqc$purity)
#'
#' subset_by_segment_totalcn(x, 2)
subset_by_segment_totalcn = function(x, totalcn)
{
  if (!inherits(x, "cnaqc"))
    stop("Not a CNAqc object in input.")

  cna_calls = x %>% CNA %>%
    dplyr::mutate(total_cn = Major + minor) %>%
    dplyr::filter(total_cn %in% totalcn)

  if (nrow(cna_calls) == 0)
    stop("There are no calls with these total copy states, cannot subset.")

  return(
    init(
      mutations = x$mutations,
      cna = cna_calls,
      purity = x$purity,
      ref = x$reference_genome,
      sample = x$sample
    )
  )
}


#' Retain clonal segments with minimum number of mutations.
#'
#' @description It retains only the set of clonal segments that have mapped
#' at least `n` mutations mapped.
#'
#' @param x A CNAqc object.
#' @param n The minimum number of mutations per segment.
#'
#' @return A CNAqc object that should have less segments.
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(mutations = example_dataset_CNAqc$mutations, cna = example_dataset_CNAqc$cna, purity = example_dataset_CNAqc$purity)
#'
#' subset_by_segment_minmutations(x, 0)
#' subset_by_segment_minmutations(x, 10)
subset_by_segment_minmutations = function(x, n = 50)
{
  if (!inherits(x, "cnaqc"))
    stop("Not a CNAqc object in input.")

  cna_calls = x %>% CNA(type = "clonal")

  removed_calls = cna_calls %>% dplyr::filter(n < !!n)
  retained_calls = cna_calls %>% dplyr::filter(n >= !!n)

  if(nrow(retained_calls) == 0)
    stop("Too strict filter -- no segments would remain!")

  if (nrow(removed_calls) > 0)
  {
    nremoved_calls = removed_calls %>% nrow()
    lmuts = removed_calls$n %>% sum
    ptmuts = 100*(lmuts/x$n_mutations)
    if(ptmuts < 1) ptmuts = '< 1'

    cli::cli_alert_info(
      "Dropping {.field {nremoved_calls}} clonal segments with less than {.field {n}} \\
      mutations - losing {.field {lmuts}}/{.field {x$n_mutations}} mutations \\
      ({.field {ptmuts}%})"
    )

    return(
      init(
        mutations = x %>% Mutations(),
        cna = retained_calls,
        purity = x$purity,
        ref = x$reference_genome,
        sample = x$sample
      )
    )
  }

  cli::cli_alert_warning("This filter has no effect on the input data.")

  return(x)

}

#' Subset mutations by minimum CCF
#'
#' @description
#'
#' Retains only a subset of mutations, if they have CCF above a cutoff. Note that
#' subclonal CNAs are lost upon application of this function.
#'
#' @param x A CNAqc object.
#' @param min_target_CCF The minimum CCF do be enforced.
#'
#' @param x A new CNAqc object with subset data.
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(mutations = example_dataset_CNAqc$mutations, cna = example_dataset_CNAqc$cna, purity = example_dataset_CNAqc$purity)
#'
#' # Original data
#' plot_data_histogram(x)
#'
#' # Original data if CCF is above 10%
#' plot_data_histogram(subset_by_minimum_CCF(x))
subset_by_minimum_CCF = function(x, min_target_CCF = 0.1)
{
  used_karyotypes = x$mutations$karyotype %>% unique

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

  observed_minima = x %>% Mutations() %>%
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
  subset_data = x$mutations %>%
    dplyr::mutate(remove = VAF <  cutoffs[karyotype]) %>%
    dplyr::filter(!remove)

  init(
    mutations = subset_data,
    cna = x$cna %>% dplyr::select(-segment_id,-n,-CCF),
    purity = x$purity,
    ref = x$reference_genome,
    sample = x$sample
  )
}

#' Subset only SNVs.
#'
#' @description This function removes all non-SNVs mutations, re-creating a new CNAqc
#' object with just SNVs. All analyses are lost.
#'
#' @param x A CNAqc object.
#' @param ref_nucleotides What reference alleles to use, default \code{c("A", "C", "T", "G")}.
#' @param alt_nucleotides What alternative alleles to use, default \code{c("A", "C", "T", "G")}.
#'
#' @return A new CNAqc dataset created with \code{init}.
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(mutations = example_dataset_CNAqc$mutations, cna = example_dataset_CNAqc$cna, purity = example_dataset_CNAqc$purity)
#'
#' subset_snvs(x)
subset_snvs = function(x,
                       ref_nucleotides = c("A", "C", "T", "G"),
                       alt_nucleotides = c("A", "C", "T", "G"))
{
  subset_mutations = x %>% Mutations %>%
    dplyr::filter(ref %in% ref_nucleotides, alt %in% alt_nucleotides)

  removed_mutations = x %>% Mutations %>%
    dplyr::filter(!(ref %in% ref_nucleotides) |
                    !(alt %in% alt_nucleotides))

  if (nrow(removed_mutations) > 0)
  {
    cli::cli_alert_success(
      "Retained {.field {nrow(subset_mutations)}} SNVs, removed {.field {nrow(removed_mutations)}}."
    )

    x = init(
      mutations = subset_mutations,
      cna = x$cna %>% dplyr::select(-segment_id,-n,-CCF),
      purity = x$purity,
      ref = x$reference_genome,
      sample = x$sample
    )
  }

  return(x)
}

# load('/Volumes/Data/Dropbox/160916_HMFregCPCT_FR10302749_FR12251983_CPCT02020357_cnaqc_object.RData')
# plot_data_histogram(x)
# w = init(x$mutations, x$cna, x$purity, x$reference_genome)
# wsb = subset_by_minimum_CCF(w, 0.1)
# wsbsnv = subset_mutations(wsb)
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
#' with data splits by chromosomes. In this way it is easy to run QC steps
#' per chromosome.
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
#' # Split of tje PCAWG object, for instance
#' split_by_chromosome(example_PCAWG)
split_by_chromosome = function(x,
                               chromosomes = paste0('chr', c(1:22, 'X', 'Y'))
                               )
{
  stopifnot(inherits(x, 'cnaqc'))

  objs = NULL
  nm = NULL

  for(chr in chromosomes)
  {
    clonal_mutations = x$mutations %>% filter(chr == !!chr)

    if((clonal_mutations %>% nrow()) == 0) next

    cli::cli_h3(chr)
    cat("\n")

    subclonal_mutations = x$cna_subclonal %>% filter(chr == !!chr) %>% pull(mutations) %>% Reduce(f = bind_rows)
    cnas = x$cna %>% filter(chr == !!chr)

    cnaqc_obj = init(
      clonal_mutations %>% bind_rows(subclonal_mutations),
      cna = cnas,
      purity = x$purity,
      ref = x$reference_genome,
      sample = x$sample
    )

    objs = append(objs, list(cnaqc_obj))
    nm = c(nm, chr)
  }

  names(objs) = nm

  return(objs)
}

