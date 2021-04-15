fortify_CNA_segments = function(x)
{
  C = colnames(x)

  required_columns = c('chr', 'from', 'to', 'Major', 'minor')

  if (!all(required_columns %in% C))
    stop("Bad Copy Number format, see the manual.")

  if (!('CCF' %in% C))
  {
    x$CCF = 1

    cli::cli_alert_warning(
      "Missing CCF column from CNA calls, adding CCF = 1 assuming clonal CNA calls."
    )
  }
  else
    {
      if(x$CCF %>% is.na() %>% any())
      {
        cli::cli_alert_warning(
          "NAs in the CCF column from CNA calls, removing those."
        )

        x = x %>% dplyr::filter(!is.na(CCF))
      }

      if((x$CCF < 1) %>% any())
      {
        cli::cli_alert_warning(
          "Subclonal CNAs are not yet supported, those segments will be removed."
        )

        x = x %>% dplyr::filter(CCF == 1) %>% distinct(chr, from, to, Major, minor, .keep_all = TRUE)
      }
    }

  x$from = enforce_numeric(x$from)
  x$to = enforce_numeric(x$to)
  x$Major = enforce_numeric(x$Major)
  x$minor = enforce_numeric(x$minor)

  if (!('length' %in% C))
  {
    x = x %>% dplyr::mutate(length = to - from)

    cli::cli_alert_warning(
      "Missing segments length from CNA calls, adding it to CNA calls."
    )
  }


  if('is_driver' %in% C)
  {
    if(!("gene" %in% C))
    {
      x = x %>% dplyr::mutate(gene = paste(chr, from, to, ref, alt, sep = ':'))
      cli::cli_alert_info("Driver annotation is present in mutation data (is_driver) but a gene column is missing, genome coordinates will be used in plotting.")
    }
    else
      cli::cli_alert_info("Driver annotation is present in mutation data (is_driver), will be used in plotting.")
  }

  # Check chromosome format reference
  if(!all(grepl("chr", x$chr)))
  {
    cli::cli_alert_warning("CNA chromosomes should be in the format 'chr*', I will add a 'chr' prefix.")
    x = x %>%
      rowwise() %>%
      mutate(chr = ifelse(
        !grepl("chr", chr),
        paste0('chr', chr),
        chr
      ))
  }


  return(tibble::as_tibble(x))
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

  if(max(c(nchar(x$alt), nchar(x$ref))) > 1)
  {
    warning(
      "You are using indels mutation data, just beware that indels count-values aree less reliable than SNVs ones ...."
    )
  }

  x$from = enforce_numeric(x$from)
  x$to = enforce_numeric(x$to)
  x$DP = enforce_numeric(x$DP)
  x$NV = enforce_numeric(x$NV)
  x$VAF = enforce_numeric(x$VAF)

  # Complement required driver information if missing
  if('is_driver' %in% colnames(x))
  {
    if(!('gene' %in% colnames(x)))
    {
      cli::cli_alert_info("Drivers are annotated, but 'gene' column is missing, using mutation location.")
      x = x %>%
        dplyr::mutate(gene = paste(chr, from, to, ref, alt, sep = ':'))
    }
  }

  # Check chromosome format reference
  if(!all(grepl("chr", x$chr)))
  {
    cli::cli_alert_warning("Mutation chromosomes should be in the format 'chr*', I will add a 'chr' prefix.")
    x = x %>%
      rowwise() %>%
      mutate(chr = ifelse(
        !grepl("chr", chr),
        paste0('chr', chr),
        chr
      ))
  }

  return(tibble::as_tibble(x))
}

enforce_numeric = function(x) {
  if (!all(is.numeric(x)))
  {
    warning("[CNAqc] Enforcing numeric for values: ", paste(head(x), collapse = ', '), ', ...')
    x = as.numeric(paste(x))
  }

  x
}

