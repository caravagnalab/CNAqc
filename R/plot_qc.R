#' Title
#'
#' @param x
#' @param digits
#'
#' @return
#' @export
#'
#' @examples
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
    CNAqc:::my_ggplot_theme(cex = .7) +
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
