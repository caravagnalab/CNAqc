get_prevalent_karyotype = function(x, rank = 1)
{
  rank_ploidy = x$cna %>%
    dplyr::group_by(minor, Major) %>%
    dplyr::summarise(n = sum(length)) %>%
    dplyr::arrange(desc(n))

  if(rank > nrow(rank_ploidy)) rank = nrow(rank_ploidy)

  return(
    paste0(
      rank_ploidy$Major[rank], ':', rank_ploidy$minor[rank]
    )
  )
}
