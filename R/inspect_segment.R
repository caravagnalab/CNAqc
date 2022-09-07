#' Plot VAFs across chromosomes.
#'
#' @description Plot VAFs across chromosomes for mutations mapping to
#' clonal simple CNAs. FIlters can be used to subset the data. Every segment has
#' a colour, and  maximum 74 segements can be plotted. The plots are
#' split by chromosome and karyotype.
#'
#' @param x A CNAqc object.
#' @param chrs The chromosome ids to use (e.g., "chr2"). All are used by default.
#' @param n Disregard chromosomes with less than `n` mutations mapped.
#' @param l Disregard chromosomes that span less than `l` nucleotide.
#'
#' @return A `ggplot2` object.
#'
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(mutations = example_dataset_CNAqc$mutations, cna = example_dataset_CNAqc$cna, purity = example_dataset_CNAqc$purity)
#'
#' # Deafault segments -- all chromsomes, at least 1KB, at least 200 mapped mutations.
#' inspect_segment(x)
#'
#' # Same as above, but only for chromosome 2
#' inspect_segment(x, chrs = 'chr2', n = 50)
inspect_segment = function(x,
                           chrs = c(paste0("chr", 1:22), 'chrX', 'chrY'),
                           n = 200,
                           l = 1000)
{
  stopifnot(inherits(x, 'cnaqc'))

  segment_ids = x %>%
    CNA() %>%
    dplyr::mutate(size = to - from) %>%
    dplyr::filter(n > !!n, size > l, chr %in% chrs) %>%
    dplyr::pull(segment_id)

  k <- length(segment_ids)

  if (k == 0)
  {
    message("No segments matching input criteria!")
    return(eplot())
  }

  hplot = x %>%
    Mutations() %>%
    # dplyr::ungroup %>%
    dplyr::filter(segment_id %in% segment_ids) %>%
    ggplot(aes(VAF)) +
    ggplot2::facet_grid(karyotype ~ chr, scales = 'free') +
    CNAqc:::my_ggplot_theme() +
    ggplot2::scale_x_continuous(breaks = c(0,  1), limits = c(0, 1)) +
    ggplot2::guides(fill = FALSE) +
    ggplot2::labs(
      title = paste0(k, ' segments, ', length(chrs), ' chromosomes'),
      subtitle = paste0('>', n, " mutations per segment, segment length >", l, ' bases.')
    )

  if (k < 74)
  {
    qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(
      RColorBrewer::brewer.pal,
      qual_col_pals$maxcolors,
      rownames(qual_col_pals)
    ))

    hplot = hplot +
      ggplot2::geom_histogram(binwidth = 0.01, ggplot2::aes(fill = segment_id)) +
      ggplot2::scale_fill_manual(values = col_vector)
  }
  else
    hplot = hplot +
    ggplot2::geom_histogram(binwidth = 0.01) +
    ggplot2::scale_fill_manual(values = col_vector)


  hplot
  # x$mutations %>%
  #   ungroup %>%
  #   filter(segment_id %in% segment_ids) %>%
  #   ggplot(aes(VAF, fill = segment_id)) +
  #   geom_histogram(binwidth = 0.01) +
  #   facet_grid(karyotype ~ chr, scales = 'free') +
  #   CNAqc:::my_ggplot_theme() +
  #   scale_x_continuous(breaks = c(0,  1), limits = c(0, 1)) +
  #   scale_fill_manual(values = col_vector) +
  #   guides(fill = FALSE)
}
