prepare_input_data = function(mutations, cna, tumour_purity)
{
  # Input data types and fortification
  stopifnot(is_tibble(mutations) | is.data.frame(mutations))
  stopifnot(is_tibble(cna) | is.data.frame(cna))
  stopifnot(tumour_purity > 0 |
              tumour_purity <= 1 | !is.na(tumour_purity))

  # Input mutations
  ref_nucleotides = c("A", "C", "T", "G")
  alt_nucleotides = c("A", "C", "T", "G")

  mutations = fortify_mutation_calls(mutations) %>%
    mutate(
      type = ifelse(
        (ref %in% ref_nucleotides) & (alt %in% alt_nucleotides),
        'SNV',
        'indel'
      )
    )

  all_mutations = mutations

  # Check specific mappability features for driver
  driver_mutations = NULL
  if(all(c("is_driver", "driver_label") %in% colnames(all_mutations)))
  {
    driver_mutations = all_mutations %>% filter(is_driver)

    if(nrow(driver_mutations) == 0) driver_mutations = NULL
  }

  if(!is.null(driver_mutations))
  {
    cli::cli_alert_success("Found annotated driver mutations: {.field {driver_mutations$driver_label}}.")
  }

  nsnvs = (all_mutations$type == "SNV") %>% sum()
  nindel = nrow(all_mutations) - nsnvs
  psnvs = round(nsnvs/(nsnvs + nindel) * 100)
  cli::cli_alert_success("Fortified calls for {.field {nrow(all_mutations)}} somatic mutations: {.field {nsnvs}} SNVs ({.field {psnvs}%}) and {.field {nindel}} indels.")

  # CNAs - clonal and subclonal
  cna_all = fortify_CNA_segments(cna)

  cna_clonal = cna_all$clonal %>% mutate(segment_id = paste(chr, from, to, Major, minor, CCF, sep = ':'))

  cna_subclonal = cna_all$subclonal
  if(!is.null(cna_subclonal)) cna_subclonal = cna_subclonal %>%
    mutate(segment_id = paste(chr, from, to, Major, minor, CCF, sep = ':'))

  ncnacl = cna_clonal %>% nrow()
  ncnasbcl = ifelse(is.null(cna_subclonal), 0, cna_subclonal %>% nrow())

  cli::cli_alert_success("Fortified CNAs for {.field {ncnacl + ncnasbcl}} segments: {.field {ncnacl}} clonal and {.field {ncnasbcl}} subclonal.")

  # Mapping mutations
  # cli::cli_alert_info("Mapping {.field {nrow(snvs)}} mutations")
  #   paste0("Input ",
  #   'n = ',
  #   nrow(snvs),
  #   " mutations for ",
  #   nrow(cna),
  #     " CNA segments (",
  #     ncnacl,
  #     " clonal, ",
  #     ncnasbcl,
  #     " subclonal)"
  #   )
  # )

  # Mapping mutations to clonal segments
  mutations = map_mutations_to_clonal_segments(mutations, cna_clonal)

  # Notify NA (not mapped)
  not_mappable = sum(is.na(mutations$karyotype))

  if(not_mappable > 0)
  {
    # cli::cli_alert_danger('{.field {not_mappable}} mutations cannot be mapped to clonal CNAs and will be removed.')
    mutations = mutations %>% dplyr::filter(!is.na(karyotype))
  }

  nmutations = nrow(mutations)

  # Tabular of mapping per segments (count)
  tab = mutations %>%
    group_by(segment_id) %>%
    summarise(n = n()) %>%
    arrange(desc(n))

  # Stats about mappability
  num_mappable = sum(!is.na(mutations$karyotype))
  perc_mappable = round(num_mappable / nsnvs * 100)

  # cat("\n")
  cli::cli_alert_success(
    "{.field {num_mappable}} mutations mapped to clonal CNAs."
  )

  # Mapping mutations subclonal segments
  if(!is.null(cna_subclonal))
    {
    cna_subclonal = map_mutations_to_subclonal_segments(mutations = all_mutations, cna_subclonal)

    if(nrow(cna_subclonal) == 0)
    {
      cna_subclonal = NULL
      cli::cli_alert_info("No mutations mapped to subclonal CNAs.")
    }
    else
    {
      nsubcl = sapply(cna_subclonal$mutations, nrow) %>% sum
      cli::cli_alert_success(
        "{.field {nsubcl}} mutations mapped to subclonal CNAs."
      )
    }
  }

  # Check if we mapped all the drivers
  idfy = function(x) { x %>% dplyr::mutate(id = paste(chr, from, to, ref, alt)) }

  if(!is.null(driver_mutations))
  {
    ids_clonal = mutations %>% idfy %>% pull(id)
    ids_subclonal = NULL

    if(!is.null(cna_subclonal))
      ids_subclonal = cna_subclonal$mutations %>%
        Reduce(f = dplyr::bind_rows) %>%
        idfy %>% pull(id)

    all_ids = c(ids_clonal, ids_subclonal)
    driver_ids = driver_mutations %>% idfy %>% pull(id)

    missing_drivers = which(!(driver_ids %in% all_ids), arr.ind = TRUE)

    if(length(missing_drivers) > 0)
    {
      missing = driver_mutations$driver_label[missing_drivers]
      missing = paste("Driver(s): ", paste(missing, collapse = ', '))

      cat("\n")
      cli::boxx(
        "Driver cannot be mapped - out of any segment!",
        header = missing,
        float = "center",
        col = 'white',
        border_col = 'black',
        background_col = 'indianred3'
      ) %>% cat
      cat("\n")

      # cli::boxx(
      #   "Lost driver mutation(s) during mappability (reported below) - out of any segment!",
      #   float = "center",
      #   col = 'white',
      #   background_col = 'red'
      #   ) %>% cat
      # cat("\n")
      #
      # driver_mutations[missing_drivers, ] %>%
      #   dplyr::select(
      #     chr, from, to, ref, alt, NV, DP, VAF,
      #     is_driver, driver_label
      #   ) %>%
      #   print()
    }

  }

  return(list(mutations = mutations, cna_clonal = cna_clonal, cna_subclonal = cna_subclonal, tab = tab))
}
