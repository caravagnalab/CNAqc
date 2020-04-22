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
  if(x$n_snvs < N) return(x)

  xx = x$snvs %>% dplyr::sample_n(N)

  xx_d = data.frame(stringsAsFactors = FALSE)
  if(keep_drivers & 'is_driver' %in% colnames(x$snvs))
    xx_d = x$snvs %>% dplyr::filter(is_driver)

  CNAqc::init(
    snvs = dplyr::bind_rows(xx_d, xx),
    cna = x$cna,
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
  if(!inherits(x, "cnaqc")) stop("Not a CNAqc object in input.")

  cna_calls = x$cna %>%
    dplyr::mutate(k = paste0(Major, ':', minor)) %>%
    dplyr::filter(k %in% karyotypes)

  if(nrow(cna_calls) == 0) stop("There are no calls with these karyotypes, cannot subset.")

  return(
    CNAqc::init(snvs = x$snvs, cna = cna_calls, purity = x$purity, ref = x$reference_genome)
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
  if(!inherits(x, "cnaqc")) stop("Not a CNAqc object in input.")

  cna_calls = x$cna %>%
    dplyr::mutate(total_cn = Major + minor) %>%
    dplyr::filter(total_cn %in% totalcn)

  if(nrow(cna_calls) == 0) stop("There are no calls with these total copy states, cannot subset.")

  return(
    CNAqc::init(snvs = x$snvs, cna = cna_calls, purity = x$purity, ref = x$reference_genome)
  )
}

