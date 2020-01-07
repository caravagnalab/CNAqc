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
    # pio::pioHdr('CNAqc - Copy Number Alteration quality check')
    # pio::pioStr("Author : ", "Giulio Caravagna <gcaravagn@gmail.com>", suffix = '\n')
    # pio::pioStr("GitHub : ", "caravagn/CNAqc", suffix = '\n')
    # pio::pioStr("   WWW : ", "https://caravagn.github.io/CNAqc/", suffix = '\n')
    #
    #
    # cat(
    #   "\n > CNAqc is part of the", crayon::green("\"evoverse\""),
    #   crayon::blue("[https://bit.ly/2orn94e]"),
    #   "- a collection of packages to implement Cancer Evolution analyses from cancer sequencing data.\n"
    # )


    pk = 'CNAqc'
    pk_l = 'Copy Number Alteration quality check'
    www = "https://caravagn.github.io/CNAqc/"
    em = "gcaravagn@gmail.com"

    cli::cli_alert_success(
      'Loading {.field {pk}}, {.emph \'{pk_l}\'}. Support : {.url { www}}' )


    options(CNAqc_welcome_message = FALSE)
  }

  invisible()
}
