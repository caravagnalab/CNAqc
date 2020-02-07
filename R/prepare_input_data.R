prepare_input_data = function(snvs, cna, tumour_purity)
{
  # Input data types and fortification
  stopifnot(is_tibble(snvs) | is.data.frame(snvs))
  stopifnot(is_tibble(cna) | is.data.frame(cna))
  stopifnot(tumour_purity > 0 |
              tumour_purity <= 1 | !is.na(tumour_purity))

  pio::pioHdr("CNAqc - CNA Quality Check")
  cat('\n')

  snvs = fortify_mutation_calls(snvs)
  cna = fortify_CNA_segments(cna) %>%
    mutate(
      segment_id = paste(chr, from, to, Major, minor, CCF, sep = ':')
      # segment_id = paste0('__sgm_', row_number())
      )

  nsnvs = nrow(snvs)
  ncna = nrow(cna)

  ncnacl = sum(cna$CCF == 1)
  ncnasbcl = sum(cna$CCF < 1)

  # Mapping mutations
  cli::cli_alert_info(
    paste0("Input ",
    'n = ',
    nsnvs,
    " mutations for ",
    ncna,
      " CNA segments (",
      ncnacl,
      " clonal, ",
      ncnasbcl,
      " subclonal)"
    )
  )


  snvs = map_mutations_to_segments(snvs, cna %>% filter(CCF == 1))

  # Tabular of mapping per segments (count)
  tab = snvs %>%
    group_by(segment_id) %>%
    summarise(n = n()) %>%
    arrange(desc(n))

  # Stats about mappability
  num_mappable = sum(!is.na(snvs$karyotype))
  perc_mappable = round(num_mappable / nsnvs * 100)

  cat("\n")
  cli::cli_alert_success(
    paste0("Mapped n = ",
    num_mappable,
    " mutations to clonal segments (", perc_mappable, '% of input)')
  )


  return(list(snvs = snvs, cna = cna, tab = tab))
}
