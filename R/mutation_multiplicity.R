mutmult_single_copy = function(x, karyotype)
{
  cli::cli_rule("Computing mutation multiplicity for single-copy karyotype {.field {karyotype}}")

  A = as.numeric(strsplit(karyotype, ':')[[1]][1])
  B = as.numeric(strsplit(karyotype, ':')[[1]][2])

  # Karyotype specific mutations - clonal segments
  cl_seg = x$cna %>%
    dplyr::filter(CCF == 1) %>%
    dplyr::pull(segment_id)

  snvs_k = x$snvs %>%
    dplyr::filter(karyotype == !!karyotype, segment_id %in% cl_seg) %>%
    dplyr::mutate(
      mutation_multiplicity = 1,
      CCF = ccf_adjustment_fun(VAF, B, A, x$purity, mutation_multiplicity)
    )

  return(list(mutations = snvs_k, params = NULL))
}

mutmult_two_copies = function(x, karyotype, entropy_quantile)
{
  cli::cli_rule(
    "Computing mutation multiplicity for karyotype {.field {karyotype}} ~ q = {.value {entropy_quantile}}"
  )

  A = as.numeric(strsplit(karyotype, ':')[[1]][1])
  B = as.numeric(strsplit(karyotype, ':')[[1]][2])

  # Karyotype specific mutations - clonal segments
  cl_seg = x$cna %>%
    dplyr::filter(CCF == 1) %>%
    dplyr::pull(segment_id)

  snvs_k = x$snvs %>%
    dplyr::filter(karyotype == !!karyotype, segment_id %in% cl_seg)

  # Expected VAF for 1 and 2 copies of the mutation
  #
  # Assumption: the aneuploidy state is immediately reached
  # out of a 1:1 state, and therefore we only care about
  # mutations in 1 copy (pre), and 2 copies (post).
  expectation = CNAqc:::expected_vaf_peak(A, B, x$purity) %>%
    mutate(label = ifelse(mutation_multiplicity == 1, "One copy", "Two copies"))

  med_coverage = median(snvs_k$DP, na.rm = TRUE)

  cli::cli_alert_info("Expected Binomial peak(s) for these calls (1 and 2 copies)")
  pioDisp(expectation)

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

  n = ceiling(med_coverage)

  # Bin(p1, n) and Bin(p2, n)
  d_1 = CNAqc:::binomial_density(p_1, n)
  d_2 = CNAqc:::binomial_density(p_2, n)

  # Then we obtain the Binomial quantile ranges for these
  # two distributions, which we use to consider only assingments
  # that have a minimum probability support
  rg_1 = CNAqc:::binomial_quantile_ranges(p_1, n, quantile_left = 0.01, quantile_right = 0.99)
  rg_2 = CNAqc:::binomial_quantile_ranges(p_2, n, quantile_left = 0.01, quantile_right = 0.99)

  # We want to create a mixture model: pi_1 * Bin(p1, n) + (1 - pi_1) * Bin(p2, n)
  # to model the mixture of those two Binomial distributions. To determine
  # the mixing proportions of this mixture we do some empirical trick of
  # get the number of observations between the two Binomial quantile ranges
  # that we have just computed.
  n_rg_1 = snvs_k %>% filter(VAF > rg_1[1], VAF < rg_1[2]) %>% nrow
  n_rg_2 = snvs_k %>% filter(VAF > rg_2[1], VAF < rg_2[2]) %>% nrow

  # Compute the actual mixing proportions, and re-scale the densities accordingly
  mixing = c(n_rg_1, n_rg_2) / (n_rg_1 + n_rg_2)
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

  filtered_profile_entropy = profile_entropy %>% dplyr::filter(diff > q)

  # Points determined cutting the H(x)
  lp = filtered_profile_entropy$x[1]
  rp = filtered_profile_entropy$x[nrow(filtered_profile_entropy)]

  cli::cli_alert_info("H(x)-derived cutoffs: [{.value {lp}}; {.value {rp}}]")

  # Assignemnts based on lp
  snvs_k = snvs_k %>%
    rowwise() %>%
    mutate(
      mutation_multiplicity = case_when(VAF < lp ~ '1',
                                             VAF > rp ~ '2',
                                             TRUE ~ 'NA'),
     mutation_multiplicity = as.numeric(mutation_multiplicity),
     CCF = ifelse(
       !is.na(mutation_multiplicity),
       ccf_adjustment_fun(VAF, B, A, x$purity, mutation_multiplicity),
       NA)
    ) %>%
    ungroup()

  return(list(
    mutations = snvs_k,
    params = list(
      expectation = expectation,
      joint = joint,
      profile_entropy = profile_entropy,
      filtered_profile_entropy = filtered_profile_entropy,
      entropy_quantile = entropy_quantile,
      p_12 = c(p_1, p_2),
      d_12 = c(d_1, d_2),
      rg_12 = c(rg_1, rg_2)
    )
  ))
}


plot_mutation_multiplicity_entropy = function(x, karyotype)
{
  computation = x$CCF_estimates[[karyotype]]
  snvs_k = computation$mutations

  # Mono-peak
  magnitude_plot = class_plot = entropy_plot = ggplot() + geom_blank()

  colors = c(`1` = 'forestgreen',
             `2` = 'darkorange3',
             `NA` = 'indianred2')

  colors = wesanderson::wes_palette("Royal1", n = 3, type = 'discrete')

  if(!(karyotype %in% c('1:1', '1:0')))
  {
    # Double peaks
    profile_entropy = computation$params$profile_entropy
    filtered_profile_entropy = computation$params$filtered_profile_entropy
    entropy_quantile = computation$params$entropy_quantile
    joint = computation$params$joint
    p_12 = computation$params$p_12

    lp = profile_entropy$x[1]
    rp = profile_entropy$x[nrow(profile_entropy)]

    q = quantile(profile_entropy$diff, entropy_quantile)

    # Magnitude plot
    magnitude_plot =
        ggplot(profile_entropy) +
        geom_histogram(aes(diff), bins = 100) +
        geom_vline(
          xintercept = q,
          color = 'darkred',
          linetype = 'dashed',
          size = .3
        ) +
        labs(
          title = bquote(partialdiff * 'H(x) ~ ' ~ hat(q) ~ '=' ~ .(entropy_quantile)),
          y = 'Magnitude',
          x = bquote(partialdiff * 'H(x)')
        ) +
        CNAqc:::my_ggplot_theme()

      class_plot =
        ggplot(joint) +
        geom_line(aes(x = x, y = A), color = 'darkred', linetype = 'dashed') +
        geom_vline(xintercept = p_12,
                   color = 'black',
                   linetype = 'dashed') +
        CNAqc:::my_ggplot_theme() +
        labs(title = paste0("Binomial Mixture"),
             y = 'Density',
             x = 'VAF')

      # Entropy plot
      entropy_plot = ggplot(joint) +
        geom_line(aes(x = x, y = entropy), color = 'black', size = .3) +
        geom_line(aes(x = x, y = diff), color = 'darkred', size = .3) +
        geom_line(
          data = profile_entropy,
          aes(x = x, y = entropy),
          color = 'black',
          size = .3
        ) +
        geom_vline(
          xintercept = c(lp, rp),
          color = 'steelblue',
          linetype = 'dashed'
        ) +
        CNAqc:::my_ggplot_theme() +
        labs(
          title = bquote("H(x) ~ I = [" * .(round(lp, 2)) * '; ' * .(round(rp, 2)) * ']'),
          subtitle = paste0(),
          y = 'Profile',
          x = bquote('H(x) (black) and' ~ partialdiff * 'H(x) (red)')
        )

    } # end: with multiple peaks

  # Mutation plots
  med_coverage = median(snvs_k$DP, na.rm = TRUE)

  mutation_plot = ggplot(snvs_k, aes(VAF, y = ..count.. / sum(..count..))) +
      geom_histogram(binwidth = 0.01, aes(fill = paste(mutation_multiplicity))) +
      xlim(0, 1) +
      CNAqc:::my_ggplot_theme() +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors) +
      guides(color = FALSE, fill = guide_legend('Copies')) +
      labs(
        y = paste0('Density (simulated ', med_coverage, 'x)'),
        title = paste0("Mutation multiplicity")
      )

  if (!(karyotype %in% c('1:1', '1:0')))
  {
    d_12 = computation$params$d_12
    rg_12 = computation$params$rg_12
    expectation = computation$params$expectation

    mutation_plot = mutation_plot +
      geom_line(
        data = data.frame(d_12[1:3]) %>% bind_rows(data.frame(d_12[2 * 1:3])) ,
        aes(x = x, y = mixture_y),
        inherit.aes = FALSE,
        size = .3
      ) +
      geom_vline(xintercept = rg_12,
                 linetype = 3,
               size = .4) +
    geom_vline(
      data = expectation,
      size = .5,
      linetype = 'dashed',
      color = 'indianred2',
      aes(xintercept = peak, color = label)
    )
  }

    # CCF plots
    CCF_plot = ggplot(snvs_k, aes(CCF, y = ..count.. / sum(..count..))) +
      geom_histogram(binwidth = 0.01, aes(fill = paste(mutation_multiplicity))) +
      CNAqc:::my_ggplot_theme() +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors) +
      guides(color = FALSE, fill = guide_legend('Copies')) +
      labs(y = paste0('Density'),
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

# Compute CCF values from mutation multiplicity
# m - minor allele
# M - Major allele
# p - purity
# mut.allele - mutation multiplicity
ccf_adjustment_fun = function(v, m, M, p, mut.allele =1)
{
  CN = m+M

  v * ((CN-2) * p + 2) / (mut.allele * p)
}

