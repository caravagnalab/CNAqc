map_mutations_to_clonal_segments = function(mutations, cna)
{
  if('karyotype' %in% colnames(mutations))
    warning("[CNAqc] a karyotype column is present in CNA calls, and will be overwritten")

  mutations$karyotype = NA
  mutations$segment_id = NA

  # pb <- progress_estimated(nrow(cna), min_time = 0)

  pb = progress::progress_bar$new(format = paste0(" ▣ :spin [:bar] :percent [ETA :eta] ▶ :elapsedfull"),
                                  total = nrow(cna), clear = TRUE, complete = "~", incomplete = " ",
                                  current = ">", width = 90)

  for(i in 1:nrow(cna))
  {
    # setTxtProgressBar(pb, i)
    # pb$tick()$print()
    pb$tick()


    mappable = which(
      mutations$chr == cna$chr[i] &
        mutations$from >= cna$from[i] &
        mutations$to <= cna$to[i]
        )

    mutations$karyotype[mappable] = paste0(cna$Major[i], ':', cna$minor[i])
    mutations$segment_id[mappable] = cna$segment_id[i]
  }

  mutations
}


map_mutations_to_subclonal_segments = function(mutations, cna_subclonal)
{
  # Map all input - retain only supported peaks
  supported_peaks = c("1:0", "1:1", "2:0", "2:1", "2:2")

  cna_subclonal = cna_subclonal %>%
    mutate(
      karyotype = paste(Major, minor, sep = ':'),
      karyotype_2 = paste(Major_2, minor_2, sep = ':'),
      analysed =
        (karyotype %in% supported_peaks) &
        (karyotype_2 %in% supported_peaks)
    )

  all_genome = (cna_subclonal$to - cna_subclonal$from) %>% sum()

  if(any(!cna_subclonal$analysed))
  {
    plost = cna_subclonal %>%
      dplyr::filter(!analysed) %>%
      mutate(length = to - from) %>%
      pull(length) %>%
      sum()
    mb = round(plost / 1e6)
    plost = round(plost/all_genome * 100)

    nlost = cna_subclonal %>%  dplyr::filter(!analysed) %>% nrow()
    ntot = cna_subclonal %>% nrow()

    cli::cli_alert_warning(
      "{.field {nlost}}/{.field {ntot}} ({.field {mb}} Mb, {.field {plost}%}) have subclonal CNAs that are not supported by CNAqc and will be filtered."
    )
    # cna_subclonal %>% dplyr::filter(!analysed) %>% print()
    cna_subclonal = cna_subclonal %>% dplyr::filter(analysed)
  }

  # Map mutations to the actual segments
  cna_subclonal$mutations = NULL

  # pb <- progress_estimated(nrow(cna_subclonal), min_time = 0)
  pb = progress::progress_bar$new(format = paste0(" ▣ :spin [:bar] :percent [ETA :eta] ▶ :elapsedfull"),
                                  total = nrow(cna_subclonal), clear = TRUE, complete = "~", incomplete = " ",
                                  current = ">", width = 90)

  for (i in 1:nrow(cna_subclonal))
  {
    # pb$tick()$print()
    pb$tick()

    cna_subclonal$mutations[i] = mutations %>%
      dplyr::filter(chr == cna_subclonal$chr[i],
                    from  >= cna_subclonal$from[i],
                    to <= cna_subclonal$to[i]) %>% list()
  }

  if (any(sapply(cna_subclonal$mutations, nrow) >= 1))
  {
    to_remove = cna_subclonal[sapply(cna_subclonal$mutations, nrow) >= 1,]
    cli::cli_alert_warning("{.field {nrow(to_remove)}} subclonal CNAs have <=1 mutations mapped and will be removed.")
    # to_remove %>% print()

    cna_subclonal = cna_subclonal[sapply(cna_subclonal$mutations, nrow) >
                                      1,]
  }

  cna_subclonal$n = sapply(cna_subclonal$mutations, nrow)

  return(cna_subclonal)
}
