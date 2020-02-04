# Empty plot
eplot = function()
{
  ggplot() +
    geom_text(data = data.frame(x = 0, y = 0, label = "NA"), aes(x=x,y=y,label=label) , size = 16) +
    theme_void()  +
    theme(plot.background = element_rect(color = 'gray', size = 1))
}