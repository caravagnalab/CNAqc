.onLoad <- function(libname, pkgname)
{
  # =-=-=-=-=-=-
  # Required packages will be listed here
  # =-=-=-=-=-=-
  requirements = c('tidyverse', 'pio', 'crayon', 'ggpubr', 'peakPick', 'RColorBrewer')

  suppressMessages(sapply(requirements, require, character.only = TRUE))

  # =-=-=-=-=-=-
  # Package options
  # =-=-=-=-=-=-
  options(pio.string_fg_colour = crayon::bgYellow$black)

  # =-=-=-=-=-=-
  # Header
  # =-=-=-=-=-=-

  CNAqc_welcome_message =  getOption('CNAqc_welcome_message', default = TRUE)

  if(CNAqc_welcome_message)
  {
    pio::pioHdr('CNAqc - Copy Number Alteration quality check')
    pio::pioStr("Author : ", "Giulio Caravagna <gcaravagn@gmail.com>", suffix = '\n')
    pio::pioStr("GitHub : ", "caravagn/CNAqc", suffix = '\n')

    options(CNAqc_welcome_message = FALSE)
  }

  invisible()
}
