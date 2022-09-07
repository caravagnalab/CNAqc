# Empty plot
eplot = function()
{
  ggplot2::ggplot(data = data.frame(x = 0, y = 0, label = "X"),
                  ggplot2::aes(x = x, y = y, label = label)) +
    # theme_void()  +
    CNAqc:::my_ggplot_theme() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(
        colour = "black",
        fill = NA,
        # size = 1,
        linetype = 'dashed'
      ),
      panel.background = ggplot2::element_rect(fill = 'gainsboro'),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank()
    )
}
