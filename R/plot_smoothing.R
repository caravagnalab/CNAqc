#' Plot smoothed and non-smoohted segments.
#'
#' @description
#'
#' Upon using the `smooth_segments` function, this function plots a multipanel figure
#' with the segments profile before and after smoothing.
#'
#' @seealso smooth_segments
#'
#' @param x A CNAq object.
#'
#' @return A \code{ggpubr} figure..
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(mutations = example_dataset_CNAqc$mutations, cna = example_dataset_CNAqc$cna, purity = example_dataset_CNAqc$purity)
#'
#' # Smooth
#' x = smooth_segments(x)
#'
#' # View differences
#' plot_smoothing(x)
plot_smoothing = function(x)
{
  # test it has smoothing
  if(all(is.null(x$before_smoothing))) {
    cli::cli_alert_danger("This CNAqc object does not have smoothing data associated.")
    return(ggplot() + geom_blank())
  }

  bpb = plot_segment_size_distribution(x$before_smoothing)
  Mb = ggplot_build(bpb)$layout$panel_scales_y[[1]]$range$range

  bpa = plot_segment_size_distribution(x)
  Ma = ggplot2::ggplot_build(bpa)$layout$panel_scales_y[[1]]$range$range

  Y_max_bplot = max(Mb, Ma)


  before = cowplot::plot_grid(
    plot_segments(x$before_smoothing) + ggplot2::labs(title = "Before smoothing"),
    bpb + ggplot2::ylim(0, Y_max_bplot),
    nrow = 1,
    ncol = 2,
    rel_widths = c(4, 1),
    align = 'h'
  )

  after = cowplot::plot_grid(
    plot_segments(x) + ggplot2::labs(title = "After smoothing"),
    bpa + ggplot2::ylim(0, Y_max_bplot),
    nrow = 1,
    ncol = 2,
    rel_widths = c(4, 1),
    align = 'h'
  )

  ggpubr::ggarrange(
    before,
    after,
    nrow = 2,
    ncol = 1
  )
}
