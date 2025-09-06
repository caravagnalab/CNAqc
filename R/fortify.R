fortify_CNA_segments = function(x)
{
  C = colnames(x)

  if(!is.data.frame(x)) stop("CNAs must be in dataframe format.")

  # CCF column missing
  if (!('CCF' %in% C))
  {
    x$CCF = 1
    C = colnames(x)

    cli::cli_alert_warning(
      "CNAs have no CCF, assuming clonal CNAs (CCF = 1)."
    )
  }

  # CCF column all NA
  if (x$CCF %>% is.na() %>% all)
  {
    x$CCF = 1
    C = colnames(x)

    cli::cli_alert_warning(
      "CNAs have NA in all CCF entries, assuming clonal CNAs (CCF = 1)."
    )
  }

  # Test for subclonal
  has_subclonal = (x$CCF < 1) %>% any(na.rm = TRUE)
  if(has_subclonal) cli::cli_alert_info("{.field {sum(x$CCF < 1, na.rm = TRUE)}} subclonal CNAs detected in the data.")

  # Check all other columns
  required_columns = c('chr', 'from', 'to', 'Major', 'minor', "CCF")
  if(has_subclonal) required_columns = c(required_columns, "Major_2", "minor_2")

  if (!all(required_columns %in% C))
    stop(cli::format_error("Bad CNA format; required columns are {.field {required_columns}}."))

  # Split
  clonal = x %>% dplyr::filter(CCF == 1) %>% distinct(chr, from, to, Major, minor, .keep_all = TRUE)

  if(has_subclonal)
    subclonal = x %>% dplyr::filter(CCF < 1) %>% distinct(chr, from, to, Major, minor, Major_2, minor_2, .keep_all = TRUE)
  else
    subclonal = NULL

  # Check NA cases - clonal
  with_NA_clonal = !complete.cases(clonal[, c('chr', 'from', 'to', 'Major', 'minor', "CCF"), drop = FALSE])

  if(any(with_NA_clonal))
  {
    nna = sum(with_NA_clonal)

    cli::cli_alert_danger("{.field {sum(with_NA_clonal)}} NA entrie(s) in clonal CNAs; segments will be removed.")

    clonal = clonal[!with_NA_clonal, ]

    if(nrow(clonal) == 0) stop("No more CNAs to work with?")
  }

  # Check NA cases - subclonal
  if(has_subclonal)
  {
    with_NA_subclonal = !complete.cases(subclonal[, c('chr', 'from', 'to', 'Major', 'minor', "CCF", "Major_2", "minor_2"), drop = FALSE])

    if(any(with_NA_subclonal))
    {
      cli::cli_alert_danger("{.field {sum(with_NA_subclonal)}} NA entrie(s) in subclonal CNAs; segments will be removed.")

      subclonal = subclonal[!with_NA_subclonal, ]
      if(nrow(subclonal) == 0) {
        has_subclonal = FALSE
        subclonal = NULL}
    }
  }


  # Enforce correct data types
  clonal$from = enforce_numeric(clonal$from)
  clonal$to = enforce_numeric(clonal$to)
  clonal$Major = enforce_numeric(clonal$Major)
  clonal$minor = enforce_numeric(clonal$minor)

  if(has_subclonal) {
    subclonal$Major_2 = enforce_numeric(subclonal$Major_2)
    subclonal$minor_2 = enforce_numeric(subclonal$minor_2)
  }

  # if(has_subclonal)
  # {
  #   x$Major_2 = enforce_numeric(x$Major_2)
  #   x$minor_2 = enforce_numeric(x$minor_2)
  # }

  if (!('length' %in% C))
  {
    clonal = clonal %>% dplyr::mutate(length = to - from)
    if(has_subclonal) subclonal  = subclonal  %>% dplyr::mutate(length = to - from)

    cli::cli_alert_warning(
      "Added segments length (in basepairs) to CNA segments."
    )
  }

  # if('is_driver' %in% C)
  # {
  #   if(!("gene" %in% C))
  #   {
  #     x = x %>% dplyr::mutate(gene = paste(chr, from, to, ref, alt, sep = ':'))
  #     cli::cli_alert_info("Driver annotation is present in mutation data (is_driver) but a gene column is missing, genome coordinates will be used in plotting.")
  #   }
  #   else
  #     cli::cli_alert_info("Driver annotation is present in mutation data (is_driver), will be used in plotting.")
  # }

  # Check chromosome format reference
  if(!all(grepl("chr", clonal$chr)))
  {
    cli::cli_alert_warning("Formatting clonal CNA chromosome annotation as 'chr*'.")

    clonal = clonal %>%
      rowwise() %>%
      mutate(chr = ifelse(
        !grepl("chr", chr),
        paste0('chr', chr),
        chr
      ))
  }

  if(!all(grepl("chr", subclonal$chr)))
  {
    cli::cli_alert_warning("Formatting subclonal CNA chromosome annotation as 'chr*'.")

    subclonal = subclonal %>%
      rowwise() %>%
      mutate(chr = ifelse(
        !grepl("chr", chr),
        paste0('chr', chr),
        chr
      ))
  }

  return(list(clonal = clonal, subclonal = subclonal))
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
    cli::cli_alert_warning(
      "Detected indels mutation (substitutions with >1 reference/alternative nucleotides)."
    )
  }

  if(x[, ..required_columns] %>% is.na() %>% any()){
    cli::cli_alert_danger(
      "NA values in some of the required mutation columns, these will be removed."
    )

    x = x[complete.cases(x[, ..required_columns]), ]

    if(nrow(x) == 0) stop("No more mutations to work with?")
  }

  x$from = enforce_numeric(x$from)
  x$to = enforce_numeric(x$to)
  x$DP = enforce_numeric(x$DP)
  x$NV = enforce_numeric(x$NV)
  x$VAF = enforce_numeric(x$VAF)

  # Complement required driver information if missing
  if('is_driver' %in% colnames(x))
  {
    if(!('driver_label' %in% colnames(x)))
    {
      cli::cli_alert_danger("Drivers are annotated, but 'driver_label' column is missing; these annotations are useless.")

      # cli::cli_alert_info("Drivers are annotated, but 'gene' column is missing, using mutation location.")
      # x = x %>%
      #   dplyr::mutate(driver_label = paste(chr, from, to, ref, alt, sep = ':'))
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

