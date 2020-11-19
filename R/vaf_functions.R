expected_vaf_peak = function(major, minor, tumour_purity)
{
  karyotype = paste0(major, ':', minor)
  tumour_ploidy = major + minor
  mutation_multiplicities = unique(c(1, major))
  
  peaks = sapply(
    mutation_multiplicities,
    ascat,
    tumour_purity = tumour_purity,
    tumour_ploidy = tumour_ploidy
  ) %>% unlist()

  data.frame(
    mutation_multiplicity = mutation_multiplicities,
    karyotype = karyotype,
    peak =  peaks,
    stringsAsFactors = FALSE
  ) %>% as_tibble()
}

ascat = function(mutation_multiplicity, tumour_purity, tumour_ploidy)
{
  expected_mutant_reads = mutation_multiplicity * tumour_purity
  expected_sequencing_depth = 2 * (1 - tumour_purity) + tumour_purity * tumour_ploidy
  expected_mutant_reads / expected_sequencing_depth
}

