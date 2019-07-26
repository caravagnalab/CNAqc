map_mutations_to_segments = function(snvs, cna)
{
  if('karyotype' %in% colnames(snvs))
    warning("[CNAqc] karyotype columns in cna will be overwritten")

  snvs$karyotype = NA

  pb = txtProgressBar(min = 1, max = nrow(cna), style = 3)

  for(i in 1:nrow(cna))
  {
    setTxtProgressBar(pb, i)

    mappable = which(
      snvs$chr == cna$chr[i] &
        snvs$from >= cna$from[i] &
        snvs$to <= cna$to[i]
        )

    snvs$karyotype[mappable] = paste0(cna$Major[i], ':', cna$minor[i])
  }

  close(pb)

  snvs
}
