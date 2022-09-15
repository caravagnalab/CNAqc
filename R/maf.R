has_MAF_annotations = function(x)
{
  grepl("MAF.", colnames(x$mutations)) %>% sum() > 3
}

augment_with_maf = function(x, maf)
{
  cli::cli_h1("Augmenting a CNAqc object with its MAF")

  cli::cli_h2("Input CNAqc object")

  x %>% print()

  # MAF loading
  cli::cli_h2("MAF input")

  MAF_input = maftools::read.maf(maf)

  cat("\n")
  print(MAF_input)

  # MAF conversion
  MAF_input = MAF_input@data %>% as_tibble()

  shared_colnames = intersect(
    x$mutations %>% colnames,
    MAF_input %>% colnames
  )

  cn = colnames(MAF_input)
  colnames(MAF_input)[!(cn %in% shared_colnames)] = paste0("MAF.", cn[!(cn %in% shared_colnames)])

  MAF_input = MAF_input %>%
    dplyr::mutate(
      chr = MAF.Chromosome,
      from = MAF.Start_Position,
      to = MAF.End_Position + 1, # re-scaling required
      ref = MAF.Reference_Allele,
      alt = ifelse(MAF.Tumor_Seq_Allele1 == ref, MAF.Tumor_Seq_Allele2, MAF.Tumor_Seq_Allele1)
    )


  x$mutations = x$mutations %>%
    dplyr::left_join(MAF_input, by = c("chr", "from", "to", "ref", "alt", shared_colnames))

  return(x)
}

as_maftools_obj = function(x, only_drivers = TRUE)
{
  stopifnot(inherits(x, 'cnaqc'))
  stopifnot(x %>% has_MAF_annotations)

  cli::cli_h1("Conversion to maftools")

  mutations_data = x %>% Mutations()
  colnames(mutations_data) = gsub('MAF.', '', mutations_data %>% colnames)

  if(only_drivers)
  {
    stopifnot(x %>% has_driver_data)

    mutations_data = mutations_data %>% dplyr::filter(is_driver)
    cli::cli_alert("Using {.field {nrow(mutations_data)}} driver mutations")

  } else
    cli::cli_alert("Using {.field {nrow(mutations_data)}} mutations")

  maftools::read.maf(mutations_data, verbose = FALSE)
}

as_maftools_cohort = function(x, only_drivers = TRUE)
{
  lapply(x, as_maftools_obj, only_drivers = only_drivers) %>%
    maftools::merge_mafs()
}
