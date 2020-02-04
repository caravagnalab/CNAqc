# Empty plot
eplot = function()
{
  ggplot(data = data.frame(x = 0, y = 0, label = "X"),
            aes(x = x, y = y, label = label)) +
    # theme_void()  +
    theme(panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 1,
      linetype = 'dashed'
    ),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()) 
}