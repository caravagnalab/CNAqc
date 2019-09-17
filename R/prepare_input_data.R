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
    mutate(id = paste0('__sgm_', row_number()))

  nsnvs = nrow(snvs)
  ncna = nrow(cna)

  ncnacl = sum(cna$CCF == 1)
  ncnasbcl = sum(cna$CCF < 1)

  # Mapping mutations
  pio::pioStr(
    "\nInput ",
    'n =',
    nsnvs,
    "mutations for",
    ncna,
    paste0(
      "CNA segments (",
      ncnacl,
      " clonal, ",
      ncnasbcl,
      " subclonal)\n\n"
    )
  )

  cat(crayon::green("Mapping mutations to clonal CNA.\n"))

  snvs = map_mutations_to_segments(snvs, cna %>% filter(CCF == 1))

  # Tabular of mapping per segments (count)
  tab = snvs %>%
    group_by(segment_id) %>%
    summarise(n = n()) %>%
    arrange(desc(n))

  # Stats about mappability
  num_mappable = sum(!is.na(snvs$karyotype))
  perc_mappable = round(num_mappable / nsnvs * 100)

  pio::pioStr(
    "\nMapping.",
    'n =',
    num_mappable,
    "mutations mapped to segments",
    paste0('(~', perc_mappable, '% of input)')
  )


  return(list(snvs = snvs, cna = cna, tab = tab))
}
