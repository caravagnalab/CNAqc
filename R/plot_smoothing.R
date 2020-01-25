# w = x
# x = smooth_segments(x)

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
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
  Ma = ggplot_build(bpa)$layout$panel_scales_y[[1]]$range$range

  Y_max_bplot = max(Mb, Ma)


  before = cowplot::plot_grid(
    plot_segments(x$before_smoothing) + labs(title = "Before smoothing"),
    bpb + ylim(0, Y_max_bplot),
    nrow = 1,
    ncol = 2,
    rel_widths = c(4, 1),
    align = 'h'
  )

  after = cowplot::plot_grid(
    plot_segments(x) + labs(title = "After smoothing"),
    bpa + ylim(0, Y_max_bplot),
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
