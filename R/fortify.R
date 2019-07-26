# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

fortify_CNA_segments = function(x)
{
  C = colnames(x)

  required_columns = c('chr', 'from', 'to', 'Major', 'minor')

  if (!all(required_columns %in% C))
    stop("Bad Copy Number format, see the manual.")

  if (!('CCF' %in% C))
  {
    x$CCF = 1

    message(
      "*** Missing CCF column from CNA calls, adding CCF = 1 assuming all calls to be clonal. ***"
    )
  }

  x$from = enforce_numeric(x$from)
  x$to = enforce_numeric(x$to)
  x$Major = enforce_numeric(x$Major)
  x$minor = enforce_numeric(x$minor)

  return(x %>% as_tibble)
  # Supported formats
  # alt_chr = c('chr', 'chromosome', 'Chromosome')
  # alt_from = c('from', 'pos', 'start', 'Start')
  # alt_to = c('to', 'end', 'End')
  #
  # w_chr = which(alt_chr %in% C)
  # w_from = which(alt_from %in% C)
  # w_to = which(alt_to %in% C)
  #
  # # Stop if any is missing
  # if(length(w_chr) == 0) stop('Chromosome in the wrong format; use any of ', paste(alt_chr, collapse = ', '))
  # if(length(w_from) == 0) stop('Segment start in the wrong format; use any of ', paste(alt_from, collapse = ', '))
  # if(length(w_to) == 0) stop('Segment end in the wrong format; use any of ', paste(alt_to, collapse = ', '))
}

fortify_mutation_calls = function(x)
{
  C = colnames(x)

  required_columns = c('chr', 'from', 'to', 'ref', 'alt', 'DP', 'NV', 'VAF')

  if (!all(required_columns %in% C))
    stop("Bad mutation format, see the manual.")

  x$from = enforce_numeric(x$from)
  x$to = enforce_numeric(x$to)
  x$DP = enforce_numeric(x$DP)
  x$NV = enforce_numeric(x$NV)
  x$VAF = enforce_numeric(x$VAF)

  return(x %>% as_tibble)
}

enforce_numeric = function(x) {
  if (!all(is.numeric(x)))
  {
    warning("[CNAqc] Enforcing numeric for values: ", paste(head(x), collapse = ', '), ', ...')
    x = as.numeric(paste(x))
  }

  x
}
