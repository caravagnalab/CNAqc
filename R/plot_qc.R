#' Plot a summary of QC results.
#'
#'  @description Results from \code{analyze_peaks} can be visualised with this
#'  function. Compared to individual karyotypes fits available with function \code{plot_peaks_analysis},
#'  this function reports sumary statistics for each karyotype, and the overall score.
#'
#' @param x An object of class \code{cnaqc}, where function \code{analyze_peaks} has
#' been computed.
#' @param digits Number of digits used to round scores.
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna,example_dataset_CNAqc$purity)
#'
#' x = analyze_peaks(x)
#' plot_qc(x)
plot_qc = function(x, digits = 4)
{
  stopifnot(inherits(x, "cnaqc"))

  with_peaks = all(!is.null(x$peaks_analysis))
  if(!with_peaks){
    stop("Input does not have peaks, see ?peaks_analysis to run peaks analysis.")
  }

  rounded_score = round(x$peaks_analysis$score, digits)

  x$peaks_analysis$matches %>%
    ggplot(aes(x = mutation_multiplicity, y = karyotype, fill = score)) +
    geom_tile(aes(width = .8, height = .8, color = weight), size = 3) +
    CNAqc:::my_ggplot_theme(cex = 1) +
    guides(
      fill = guide_colorbar(bquote(rho ~ ' '), barwidth = unit(2, 'cm')),
      color = guide_colorbar(bquote(omega ~ ' '), barwidth = unit(2, 'cm'))
    ) +
    scale_fill_gradient2(high = 'steelblue') +
    scale_color_distiller(palette = 'Greens', direction = 1, limits = c(0, 1)) +
    scale_alpha_continuous() +
    geom_text(
      data = x$peaks_analysis$matches,
      aes(label = round(score, digits), color = score),
      size = 3
    ) +
    geom_text(
      data = x$peaks_analysis$matches %>% mutate(cross = ifelse(matched, '', 'X')),
      aes(label = cross),
      color = 'red'
    ) +
    labs(
      x = 'Mutation multiplicuty',
      y = 'Karyotype',
      title = "Summary matches",
      subtitle = bquote("QC score " ~ widehat(rho)  ~' = ' ~ sum(omega[i] ~ rho[i], i, )  ~' = ' ~ .(rounded_score))
    )
}
