# Empty plot
eplot = function()
{
  ggplot(data = data.frame(x = 0, y = 0, label = "X"),
            aes(x = x, y = y, label = label)) +
    # theme_void()  +
    CNAqc:::my_ggplot_theme() +
    theme(panel.border = element_rect(
      colour = "black",
      fill = NA,
      # size = 1,
      linetype = 'dashed'
    ),
    panel.background = element_rect(fill = 'gainsboro'),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank())
}
