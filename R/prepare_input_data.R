prepare_input_data = function(snvs, cna, tumour_purity)
{
  # Input data types and fortification
  stopifnot(is_tibble(snvs) | is.data.frame(snvs))
  stopifnot(is_tibble(cna) | is.data.frame(cna))
  stopifnot(tumour_purity > 0 |
              tumour_purity <= 1 | !is.na(tumour_purity))

  nsnvs = nrow(snvs)
  ncna = nrow(cna)

  pio::pioHdr("CNAqc - CNA Quality Check")
  cat('\n')

  snvs = CNAqc:::fortify_mutation_calls(snvs)
  cna = CNAqc:::fortify_CNA_segments(cna)

  ncnacl = sum(cna$CCF == 1)
  ncnasbcl = sum(cna$CCF < 1)

  # Mapping mutations
  pio::pioStr("\nInput.",
              'n =',
              nsnvs,
              "mutations for",
              ncna,
              "CNA segments (", ncnacl, "clonal, ", ncnasbcl, "subclonal)\n\n")
  snvs = CNAqc:::map_mutations_to_segments(snvs, cna %>% filter(CCF == 1))

  return(list(snvs = snvs, cna =cna))
}
