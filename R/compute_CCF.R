compute_CCF = function(snvs, cna, tumour_purity, karyotypes = c('2:1', '2:0', '2:2'))
{
  input = prepare_input_data(snvs, cna, tumour_purity)
  snvs = input$snvs
  cna = input$cna

  nk = table(snvs$karyotype)[karyotypes]

  karyotypes = names(nk)[!is.na(nk)]

  fits = lapply(karyotypes, mutation_multipl, snvs = snvs)
  names(fits) = karyotypes


  # mutation_plot + xlim(rg_1)

  # library(ggforce)
  #
  # d1_ylim = c(0, max(d_2$y))
  # mutation_plot +
  #   facet_zoom(xlim = rg_2, ylim = d1_ylim, horizontal = T) +
  #   facet_zoom(xlim = rg_1, ylim = d1_ylim, horizontal = T)

}

data('example_dataset_CNAqc')
x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna, example_dataset_CNAqc$purity)

mutation_multipl = function(x, karyotype)
{
  A = as.numeric(strsplit(karyotype, ':')[[1]][1])
  B = as.numeric(strsplit(karyotype, ':')[[1]][2])

  # Karyotype specific mutations
  snvs_k = x$snvs %>%
    filter(karyotype == !!karyotype)

  # Expected VAF for 1 and 2 copies of the mutation
  #
  # Assumption: the aneuploidy state is immediately reached
  # out of a 1:1 state, and therefore we only care about
  # mutations in 1 copy (pre), and 2 copies (post).
  expectation = CNAqc:::expected_vaf_peak(A, B, x$purity) %>%
    mutate(
      label = ifelse(mutation_multiplicity == 1, "One copy", "Two copies")
    )
  med_coverage = median(snvs_k$DP, na.rm = TRUE)

  mutation_plot = ggplot(snvs_k, aes(VAF, y = ..count.. / sum(..count..))) +
    geom_histogram(binwidth = 0.01, fill = 'gainsboro') +
    xlim(0, 1) +
    geom_vline(
      data = expectation,
      size = .3,
      linetype = 'dashed',
      aes(xintercept = peak, color = label)
    ) +
    my_ggplot_theme() +
    scale_color_brewer(palette = "Set1") +
    guides(color = guide_legend("Occurrence")) +
    labs(
      y = 'Density',
      subtitle = paste0('Coverage: ', med_coverage, 'x'),
      title = paste0("Mutation multiplicity for ", karyotype))


  # Density of 2 Binomials with peaks at p_1 and p_2
  # - first one, events before aneuploidy
  # - second one, events after aneuploidy
  domain = seq(0, 1, 0.01)
  NV = round(domain * med_coverage)

  # Before - density values
  p_1 = expectation$peak[1]
  d_1 = data.frame(
    x = NV/med_coverage,
    y = dbinom(NV, size = med_coverage, prob = p_1),
    mutation_multiplicity = 1,
    label = 'One copy'
  )
  # quantile ranges
  rg_1 = c(
    qbinom(0.01, med_coverage, p_1)/med_coverage,
    qbinom(0.99, med_coverage, p_1)/med_coverage
  )

  # After - density values
  p_2 = expectation$peak[2]
  d_2 = data.frame(
    x = NV/med_coverage,
    y = dbinom(NV, size = med_coverage, prob = p_2),
    mutation_multiplicity = 2,
    label = 'Two copies'
  )
  # quantile ranges
  rg_2 = c(
    qbinom(0.01, med_coverage, p_2)/med_coverage,
    qbinom(0.99, med_coverage, p_2)/med_coverage
  )

  # Partition the data, get the number of observations between the
  # two Binomial quantile ranges that we have just computed. Use
  # the lower left-most and the upper righmost values.
  n_rg_1 = snvs_k %>% filter(VAF > rg_1[1], VAF < rg_1[2]) %>% nrow
  n_rg_2 = snvs_k %>% filter(VAF > rg_2[1], VAF < rg_2[2]) %>% nrow

  # Compute the actual mixing proportions, and re-scale the densities accordingly
  mixing = c(n_rg_1, n_rg_2)/(n_rg_1 + n_rg_2)
  d_1$y = d_1$y * mixing[1]
  d_2$y = d_2$y * mixing[2]

  # d = d_1 %>%
  #   rename(y1 = y) %>%
  #   select(x, y1) %>%
  #   bind_cols(d_2 %>%
  #               rename(y2 = y) %>%
  #               select(y2)) %>%
  #   mutate(r = y1 >= y2) %>%
  #   as_tibble()

  den = function(x)
  {
    d1 = mixing[1] * dbinom(round(x * med_coverage), size = med_coverage, prob = p_1)
    d2 = mixing[2] * dbinom(round(x * med_coverage), size = med_coverage, prob = p_2)

    return(ifelse(d1 >= d2, 1, 2))
  }

  snvs_k = snvs_k %>%
    mutate(
      mutation_multiplicity = den(VAF),
      label = ifelse(mutation_multiplicity == 2, "Two copies", "One copy")
    )

  mutation_plot = mutation_plot +
    geom_histogram(data = snvs_k, binwidth = 0.01, aes(color = label), fill = NA, size = .1) +
    geom_line(data = d_1, aes(x = x, y = y, color = label), inherit.aes = FALSE) +
    geom_line(data = d_2, aes(x = x, y = y, color = label), inherit.aes = FALSE) +
    geom_vline(xintercept = rg_1, linetype = 3, size = .2) +
    geom_vline(xintercept = rg_2, linetype = 3, size = .2)

  # zoom
  zoom_plot = ggplot(snvs_k, aes(VAF, fill = label)) +
    geom_histogram(binwidth = 0.01) +
    facet_wrap(~label, ncol = 1, scales = 'free') +
    theme_minimal() +
    theme(legend.position = 'bottom',
          legend.key.height = unit(3, 'mm')) +
    scale_fill_brewer(palette = "Set1") +
    guides(fill = FALSE) +
    labs(
      y = 'Counts'
      )

  mutation_plot = cowplot::plot_grid(
    mutation_plot,
    zoom_plot,
    align = 'h',
    ncol = 2, scale = c(1, .8)
    )

  return(list(snvs = snvs_k, plot = mutation_plot))
}


binomial_range = function(p, n)
{
  # Domain
  domain = seq(0, 1, 0.01)

  # Bernoulli succesfull trials
  NV = round(domain * n)

  # Before - density values
  p_1 = expectation$peak[1]
  d_1 = data.frame(
    x = NV/med_coverage,
    y = dbinom(NV, size = med_coverage, prob = p_1),
    mutation_multiplicity = 1,
    label = 'One copy'
  )


}
