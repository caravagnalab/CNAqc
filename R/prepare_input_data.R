prepare_input_data = function(snvs, cna, tumour_purity)
{
  # Input data types and fortification
  stopifnot(is_tibble(snvs) | is.data.frame(snvs))
  stopifnot(is_tibble(cna) | is.data.frame(cna))
  stopifnot(tumour_purity > 0 |
              tumour_purity <= 1 | !is.na(tumour_purity))

  # Input mutations
  ref_nucleotides = c("A", "C", "T", "G")
  alt_nucleotides = c("A", "C", "T", "G")
  
  snvs = fortify_mutation_calls(snvs) %>% 
    mutate(
      type = ifelse(
        (ref %in% ref_nucleotides) & (alt %in% alt_nucleotides),
        'SNV',
        'indel'
      )
    )

  cna = fortify_CNA_segments(cna) %>%
    mutate(
      segment_id = paste(chr, from, to, Major, minor, CCF, sep = ':')
      # segment_id = paste0('__sgm_', row_number())
      )

  ncnacl = sum(cna$CCF == 1)
  ncnasbcl = sum(cna$CCF < 1)

  # Mapping mutations
  cli::cli_alert_info(
    paste0("Input ",
    'n = ',
    nrow(snvs),
    " mutations for ",
    nrow(cna),
      " CNA segments (",
      ncnacl,
      " clonal, ",
      ncnasbcl,
      " subclonal)"
    )
  )

  # Subclonal CNA calls -- raise informative warning
  if(ncnasbcl > 0)
    cli::boxx("Subclonal (CCF < 1) CNA calls are in the data.
They will be stored in the 'cna' slot of the object, but not be used for the pick detection analysis ", col = 'red')

  snvs = map_mutations_to_segments(snvs, cna %>% filter(CCF == 1))

  # Notify NA (not mapped)
  not_mappable = sum(is.na(snvs$karyotype))

  if(not_mappable > 0)
  {
    cli::cli_alert_danger('n = {.value {not_mappable}} mutations  cannot be mapped to segments and will be removed.')
    snvs = snvs %>% dplyr::filter(!is.na(karyotype))
  }

  nsnvs = nrow(snvs)
  ncna = nrow(cna)

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
