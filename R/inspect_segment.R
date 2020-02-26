#' Plot the VAF of a set of segments
#'
#' @description Plot the VAF of a set of segments that map to some
#' chromosome id, with a minimum number of mapped mutations and length.
#' THe plot is a histogram of the VAF, coloured per segement. A maximum
#' of 74 segements can be plotted in colour, otherwise colors are not used
#' and the histogram is black. The plots are facetted by chromosome id and
#' karyotype.
#'
#' @param x An object of class \code{cnaqc}.
#' @param chrs The chromosome ids to use (e.g., "chr2"). All are used by default.
#' @param n Disregard chromosome with less than `n` mutations mapped.
#' @param l Disregard chromosome that span less than `l` nucleotide.
#'
#' @return A ggplot object of the VAF histogram per segment.
#'
#' @import RColorBrewer
#'
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna,example_dataset_CNAqc$purity)
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


  segment_ids = x$cna %>%
    mutate(size = to - from) %>%
    filter(n > !!n, size > l, chr %in% chrs) %>%
    pull(segment_id)

  k <- length(segment_ids)

  if(k == 0)
  {
    message("No segments matching input criteria!")
    return(ggplot() + geom_blank())
  }

  hplot = x$snvs %>%
    ungroup %>%
    filter(segment_id %in% segment_ids) %>%
    ggplot(aes(VAF)) +
    facet_grid(karyotype ~ chr, scales = 'free') +
    CNAqc:::my_ggplot_theme() +
    scale_x_continuous(breaks = c(0,  1), limits = c(0, 1)) +
    guides(fill = FALSE) +
    labs(
      title = paste0(k, ' segments, ', length(chrs), ' chromosomes'),
      subtitle = paste0('>',n, " mutations per segment, segment length >", l, ' bases.')
    )

  if(k < 74)
  {
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual', ]
    col_vector = unlist(mapply(
      brewer.pal,
      qual_col_pals$maxcolors,
      rownames(qual_col_pals)
    ))

    hplot = hplot +
      geom_histogram(binwidth = 0.01, aes(fill = segment_id)) +
      scale_fill_manual(values = col_vector)
  }
  else
    hplot = hplot +
      geom_histogram(binwidth = 0.01) +
      scale_fill_manual(values = col_vector)


  hplot
  # x$snvs %>%
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
