map_mutations_to_segments = function(snvs, cna)
{
  if('karyotype' %in% colnames(snvs))
    warning("[CNAqc] a karyotype column is present in CNA calls, and will be overwritten")

  snvs$karyotype = NA
  snvs$segment_id = NA

  pb <- progress_estimated(nrow(cna), min_time = 2)

  for(i in 1:nrow(cna))
  {
    # setTxtProgressBar(pb, i)
    pb$tick()$print()

    mappable = which(
      snvs$chr == cna$chr[i] &
        snvs$from >= cna$from[i] &
        snvs$to <= cna$to[i]
        )

    snvs$karyotype[mappable] = paste0(cna$Major[i], ':', cna$minor[i])
    snvs$segment_id[mappable] = cna$segment_id[i]
  }

  snvs
}
