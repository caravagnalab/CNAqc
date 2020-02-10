#' Randomly subsample a dataset
#'
#'
#' @param x
#' @param N
#' @param keep_drivers
#'
#' @return
#' @export
#'
#' @examples
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
