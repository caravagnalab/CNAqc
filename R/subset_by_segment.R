#' Subset calls by segments' karyotype.
#'
#' @description
#'
#' @param x
#' @param karyotype
#'
#' @return
#' @export
#'
#' @examples
subset_by_segment_karyotype = function(x, karyotype)
{
  stopifnot(inherits(x, "cnaqc"), "Not a CNAqc object in input.")

  s_karyotype = strsplit(karyotype, split = ':')[[1]]

  minor = s_karyotype[2]
  Major = s_karyotype[1]

  cna_calls = x$cna %>% dplyr::filter(minor == !!minor, Major == !!Major)

  return(
    CNAqc::init(snvs = x$snvs, cna = cna_calls, purity = x$purity, ref = x$reference_genome)
  )
}
