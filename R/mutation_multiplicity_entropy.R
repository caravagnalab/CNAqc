mutation_multiplicity_entropy = function(x, karyotype, entropy_quantile = .9)
{
  # pio::pioTit("Computing mutation multiplicity with entropy H(x)")
  cli::cli_rule("Computing mutation multiplicity for karyotype {.field {karyotype}} ~ q = {.value {entropy_quantile}}")

  A = as.numeric(strsplit(karyotype, ':')[[1]][1])
  B = as.numeric(strsplit(karyotype, ':')[[1]][2])

  colors = c(
    `1` = 'forestgreen',
    `2` = 'darkorange3',
    `NA` = 'indianred2'
  )

  # Karyotype specific mutations - clonal segments
  cl_seg = x$cna %>%
    filter(CCF == 1) %>%
    pull(segment_id)

  snvs_k = x$snvs %>%
    filter(karyotype == !!karyotype, segment_id %in% cl_seg)

  # Expected VAF for 1 and 2 copies of the mutation
  #
  # Assumption: the aneuploidy state is immediately reached
  # out of a 1:1 state, and therefore we only care about
  # mutations in 1 copy (pre), and 2 copies (post).
  expectation = CNAqc:::expected_vaf_peak(A, B, x$purity) %>%
    mutate(label = ifelse(mutation_multiplicity == 1, "One copy", "Two copies"))

  med_coverage = median(snvs_k$DP, na.rm = TRUE)

  cli::cli_alert_info("Expected Binomial peaks for these calls")
  # pio::pioStr("Binomial peaks", suffix = '\n')
  pioDisp(expectation)
  # pio::pioStr("H(x) quantile", entropy_quantile, suffix = '\n')
  # pio::pioStr(" Median depth", med_coverage, suffix = '\n')


  # =-=-=-=-=-=-=-=-=-=-=-
  # Entropy-derived heuristic for the detection of points
  # that are difficult to assign
  # =-=-=-=-=-=-=-=-=-=-=-

  # We build 2 template Binomial densities to capture:
  #
  # - Bin(p1, n), events before aneuploidy
  # - Bin(p2, n), events after aneuploidy
  #
  # In both cases we take as overal number of trials (n)
  # the median coverage, and use for the success parameters
  # p1 and p2 the expected peaks as of ASCAT equation.
  #
  # Assumptions:
  # - overdispersion is small to justify a Binomial instead
  #   of a Beta-Binomial model;
  # - trials are well-represented with the median coverage;
  p_1 = expectation$peak[1]
  p_2 = expectation$peak[2]

  n = med_coverage

  # Bin(p1, n)
  d_1 = CNAqc:::binomial_density(p_1, n)

  # Bin(p2, n)
  d_2 = CNAqc:::binomial_density(p_2, n)

  # Then we obtain the Binomial quantile ranges for these
  # two distributions, which we use to consider only assingments
  # that have a minimum probability support
  rg_1 = CNAqc:::binomial_quantile_ranges(p_1, n, quantile_left = 0.01, quantile_right = 0.99)
  rg_2 = CNAqc:::binomial_quantile_ranges(p_2, n, quantile_left = 0.01, quantile_right = 0.99)

  # We want to create a mixture model
  #
  # pi_1 * Bin(p1, n) + (1 - pi_1) * Bin(p2, n)
  #
  # to model the mixture of those two Binomial distributions. To determine
  # the mixing proportions of this mixture we do some empirical trick of
  # get the number of observations between the two Binomial quantile ranges
  # that we have just computed.
  n_rg_1 = snvs_k %>% filter(VAF > rg_1[1], VAF < rg_1[2]) %>% nrow
  n_rg_2 = snvs_k %>% filter(VAF > rg_2[1], VAF < rg_2[2]) %>% nrow

  # Compute the actual mixing proportions, and re-scale the densities accordingly
  mixing = c(n_rg_1, n_rg_2)/(n_rg_1 + n_rg_2)
  d_1$mixture_y = d_1$y * mixing[1]
  d_2$mixture_y = d_2$y * mixing[2]

  cli::cli_alert_info("Mixing pre/ post aneuploidy: {.value {round(mixing, 2)}}")

  # Now we need to decide how to assign a point in order to determine the actual mutation
  # multeplicity. We want this to be using the entropy of a 2-class model, and the
  # magnitude of the differential of the entropy
  joint = CNAqc:::entropy_profile_2_class(d_1, d_2)

  # q-Quantile of the entropy that we try to cut out, profiled in between the
  # two Binomial peaks (fluctuations outside are generally irrelevant)
  profile_entropy = joint %>%
    filter(x > p_1, x < p_2)

  q = quantile(profile_entropy$diff, entropy_quantile)

  # Magnitude plot
  magnitude_plot =
    ggplot(profile_entropy) +
    geom_histogram(aes(diff), bins = 100) +
    geom_vline(xintercept = q, color = 'darkred', linetype = 'dashed', size = .3) +
    labs(title = bquote(partialdiff * 'H(x) ~ ' ~ hat(q) ~ '=' ~ .(entropy_quantile)),
         y = 'Magnitude',
         x = bquote(partialdiff * 'H(x)')) +
    my_ggplot_theme()

  profile_entropy = profile_entropy %>% filter(diff > q)

  # Points determined cutting the H(x)
  lp = profile_entropy$x[1]
  rp = profile_entropy$x[nrow(profile_entropy)]

  cli::cli_alert_info("H(x)-derived cutoffs: [{.value {lp}}; {.value {rp}}]")

  y_lp = joint$entropy[round(lp, 2) * 100]
  y_rp = joint$entropy[round(rp, 2) * 100]

  class_plot =
    ggplot(joint) +
    geom_line(aes(x = x, y = A), color = 'darkred', linetype = 'dashed') +
    geom_vline(xintercept = c(p_1, p_2), color = 'black', linetype = 'dashed') +
    my_ggplot_theme() +
    labs(title = paste0("Binomial Mixture"),
         y = 'Density',
         x = 'VAF')

  # Entropy plot
  entropy_plot = ggplot(joint) +
    geom_line(aes(x = x, y = entropy), color = 'black', size = .3) +
      geom_line(aes(x = x, y = diff), color = 'darkred', size = .3) +
      geom_line(data = profile_entropy, aes(x = x, y = entropy), color = 'black', size = .3) +
    geom_vline(xintercept = c(lp, rp), color = 'steelblue', linetype = 'dashed') +
    my_ggplot_theme() +
    labs(title = bquote("H(x) ~ I = [" * .(round(lp, 2)) * '; ' * .(round(rp, 2)) * ']' ),
         subtitle = paste0(),
         y = 'Profile',
         x = bquote('H(x) (black) and' ~ partialdiff *'H(x) (red)'))

  # Assignemnts based on lp
  snvs_k = snvs_k %>%
    mutate(
      mutation_multiplicity = case_when(
        VAF < lp ~ '1',
        VAF > rp ~ '2',
        TRUE ~ 'NA')
    )
  snvs_k$mutation_multiplicity = as.numeric(snvs_k$mutation_multiplicity)

  # Compute CCF values from mutation multiplicity
  # m - minor allele
  # M - Major allele
  # p - purity
  # mut.allele - mutation multiplicity
  adjustment_fun = function(v, m, M, p, mut.allele =1)
  {
    CN = m+M

    v * ((CN-2) * p + 2) / (mut.allele * p)
  }

  snvs_k = snvs_k %>%
    mutate(
      CCF = adjustment_fun(VAF, B, A, x$purity, mutation_multiplicity)
    )

  # Mutation plots
  mutation_plot = ggplot(snvs_k, aes(VAF, y = ..count.. / sum(..count..))) +
    geom_histogram(binwidth = 0.01, aes(fill = paste(mutation_multiplicity))) +
    xlim(0, 1) +
    geom_vline(
      data = expectation,
      size = .3,
      linetype = 'dashed',
      aes(xintercept = peak, color = label)
    ) +
    my_ggplot_theme() +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    guides(color = FALSE, fill = guide_legend('Copies')) +
    labs(
      y = paste0('Density (n = ', med_coverage, 'x)'),
      title = paste0("Mutation multiplicity"))

  mutation_plot = mutation_plot +
    geom_histogram(data = snvs_k, binwidth = 0.01, fill = NA, size = .1) +
    geom_line(data = d_1, aes(x = x, y = mixture_y), inherit.aes = FALSE, size = .3) +
    geom_line(data = d_2, aes(x = x, y = mixture_y), inherit.aes = FALSE, size = .3) +
    geom_vline(xintercept = rg_1, linetype = 3, size = .2) +
    geom_vline(xintercept = rg_2, linetype = 3, size = .2)

  # CCF plots
  CCF_plot = ggplot(snvs_k, aes(CCF, y = ..count.. / sum(..count..))) +
    geom_histogram(binwidth = 0.01, aes(fill = paste(mutation_multiplicity))) +
    my_ggplot_theme() +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    guides(color = FALSE, fill = guide_legend('Copies')) +
    labs(
      y = paste0('Density'),
      title = paste0("CCF for ", karyotype))

  figure = cowplot::plot_grid(
    CCF_plot,
    mutation_plot,
    class_plot,
    entropy_plot,
    magnitude_plot,
    # cowplot::plot_grid(magnitude_plot, entropy_plot, align = 'v', nrow = 2),
    align = 'h',
    axis = 'b',
    ncol = 5
    # scale = c(1, .8)
  )

  return(list(mutations = snvs_k, plot = figure))
}
