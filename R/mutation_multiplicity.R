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

# OLD IMPLEMENTATION -- removed with the automatic peak detection
# mutmult_two_copies = function(x, karyotype, entropy_quantile)
# {
#   cli::cli_rule(
#     "Computing mutation multiplicity for karyotype {.field {karyotype}} ~ q = {.value {entropy_quantile}}"
#   )
#
#   A = as.numeric(strsplit(karyotype, ':')[[1]][1])
#   B = as.numeric(strsplit(karyotype, ':')[[1]][2])
#
#   # Karyotype specific mutations - clonal segments
#   cl_seg = x$cna %>%
#     dplyr::filter(CCF == 1) %>%
#     dplyr::pull(segment_id)
#
#   snvs_k = x$snvs %>%
#     dplyr::filter(karyotype == !!karyotype, segment_id %in% cl_seg)
#
#   # Expected VAF for 1 and 2 copies of the mutation
#   #
#   # Assumption: the aneuploidy state is immediately reached
#   # out of a 1:1 state, and therefore we only care about
#   # mutations in 1 copy (pre), and 2 copies (post).
#   expectation = CNAqc:::expected_vaf_peak(A, B, x$purity) %>%
#     mutate(label = ifelse(mutation_multiplicity == 1, "One copy", "Two copies"))
#
#   med_coverage = median(snvs_k$DP, na.rm = TRUE)
#
#   cli::cli_alert_info("Expected Binomial peak(s) for these calls (1 and 2 copies)")
#   pioDisp(expectation)
#
#   # =-=-=-=-=-=-=-=-=-=-=-
#   # Entropy-derived heuristic for the detection of points
#   # that are difficult to assign
#   # =-=-=-=-=-=-=-=-=-=-=-
#   # We build 2 template Binomial densities to capture:
#   #
#   # - Bin(p1, n), events before aneuploidy
#   # - Bin(p2, n), events after aneuploidy
#   #
#   # In both cases we take as overal number of trials (n)
#   # the median coverage, and use for the success parameters
#   # p1 and p2 the expected peaks as of ASCAT equation.
#   #
#   # Assumptions:
#   # - overdispersion is small to justify a Binomial instead
#   #   of a Beta-Binomial model;
#   # - trials are well-represented with the median coverage;
#   p_1 = expectation$peak[1]
#   p_2 = expectation$peak[2]
#
#   n = ceiling(med_coverage)
#
#   # Bin(p1, n) and Bin(p2, n)
#   d_1 = CNAqc:::binomial_density(p_1, n, N_bins = 1000)
#   d_2 = CNAqc:::binomial_density(p_2, n, N_bins = 1000)
#
#   # Then we obtain the Binomial quantile ranges for these
#   # two distributions, which we use to consider only assingments
#   # that have a minimum probability support
#   rg_1 = CNAqc:::binomial_quantile_ranges(p_1, n, quantile_left = 0.01, quantile_right = 0.99)
#   rg_2 = CNAqc:::binomial_quantile_ranges(p_2, n, quantile_left = 0.01, quantile_right = 0.99)
#
#   # We want to create a mixture model: pi_1 * Bin(p1, n) + (1 - pi_1) * Bin(p2, n)
#   # to model the mixture of those two Binomial distributions. To determine
#   # the mixing proportions of this mixture we do some empirical trick of
#   # get the number of observations between the two Binomial quantile ranges
#   # that we have just computed.
#   n_rg_1 = snvs_k %>% filter(VAF > rg_1[1], VAF < rg_1[2]) %>% nrow
#   n_rg_2 = snvs_k %>% filter(VAF > rg_2[1], VAF < rg_2[2]) %>% nrow
#
#   # Compute the actual mixing proportions, and re-scale the densities accordingly
#   mixing = c(n_rg_1, n_rg_2) / (n_rg_1 + n_rg_2)
#   d_1$mixture_y = d_1$y * mixing[1]
#   d_2$mixture_y = d_2$y * mixing[2]
#
#   cli::cli_alert_info("Mixing pre/ post aneuploidy: {.value {round(mixing, 2)}}")
#
#   # Now we need to decide how to assign a point in order to determine the actual mutation
#   # multeplicity. We want this to be using the entropy of a 2-class model, and the
#   # magnitude of the differential of the entropy
#   joint = CNAqc:::entropy_profile_2_class(d_1, d_2)
#
#   # q-Quantile of the entropy that we try to cut out, profiled in between the
#   # two Binomial peaks (fluctuations outside are generally irrelevant)
#   profile_entropy = joint %>%
#     filter(x > p_1, x < p_2)
#
#   q = quantile(profile_entropy$diff, entropy_quantile)
#
#   filtered_profile_entropy = profile_entropy %>% dplyr::filter(diff > q)
#
#   # Points determined cutting the H(x)
#   lp = filtered_profile_entropy$x[1]
#   rp = filtered_profile_entropy$x[nrow(filtered_profile_entropy)]
#
#   cli::cli_alert_info("H(x)-derived cutoffs: [{.value {lp}}; {.value {rp}}]")
#
#   # Assignemnts based on lp
#   snvs_k = snvs_k %>%
#     rowwise() %>%
#     mutate(
#       mutation_multiplicity = case_when(VAF < lp ~ '1',
#                                              VAF > rp ~ '2',
#                                              TRUE ~ 'NA'),
#      mutation_multiplicity = as.numeric(mutation_multiplicity),
#      CCF = ifelse(
#        !is.na(mutation_multiplicity),
#        ccf_adjustment_fun(VAF, B, A, x$purity, mutation_multiplicity),
#        NA)
#     ) %>%
#     ungroup()
#
#   return(list(
#     mutations = snvs_k,
#     params = list(
#       expectation = expectation,
#       joint = joint,
#       profile_entropy = profile_entropy,
#       filtered_profile_entropy = filtered_profile_entropy,
#       entropy_quantile = entropy_quantile,
#       p_12 = c(p_1, p_2),
#       d_12 = c(d_1, d_2),
#       rg_12 = c(rg_1, rg_2)
#     )
#   ))
# }
#
# plot_mutation_multiplicity_entropy = function(x, karyotype)
# {
#   computation = x$CCF_estimates[[karyotype]]
#   snvs_k = computation$mutations
#
#   # Mono-peak
#   magnitude_plot = class_plot = entropy_plot = ggplot() + geom_blank()
#
#   colors = c(`1` = 'forestgreen',
#              `2` = 'darkorange3',
#              `NA` = 'indianred2')
#
#   colors = wesanderson::wes_palette("Royal1", n = 3, type = 'discrete')
#
#   if(!(karyotype %in% c('1:1', '1:0')))
#   {
#     # Double peaks
#     profile_entropy = computation$params$profile_entropy
#     filtered_profile_entropy = computation$params$filtered_profile_entropy
#     entropy_quantile = computation$params$entropy_quantile
#     joint = computation$params$joint
#     p_12 = computation$params$p_12
#
#     lp = profile_entropy$x[1]
#     rp = profile_entropy$x[nrow(profile_entropy)]
#
#     q = quantile(profile_entropy$diff, entropy_quantile)
#
#     # Magnitude plot
#     magnitude_plot =
#         ggplot(profile_entropy) +
#         geom_histogram(aes(diff), bins = 100) +
#         geom_vline(
#           xintercept = q,
#           color = 'darkred',
#           linetype = 'dashed',
#           size = .3
#         ) +
#         labs(
#           title = bquote(partialdiff * 'H(x) ~ ' ~ hat(q) ~ '=' ~ .(entropy_quantile)),
#           y = 'Magnitude',
#           x = bquote(partialdiff * 'H(x)')
#         ) +
#         CNAqc:::my_ggplot_theme()
#
#       class_plot =
#         ggplot(joint) +
#         geom_line(aes(x = x, y = A), color = 'darkred', linetype = 'dashed') +
#         geom_vline(xintercept = p_12,
#                    color = 'black',
#                    linetype = 'dashed') +
#         CNAqc:::my_ggplot_theme() +
#         labs(title = paste0("Binomial Mixture"),
#              y = 'Density',
#              x = 'VAF')
#
#       # Entropy plot
#       entropy_plot = ggplot(joint) +
#         geom_line(aes(x = x, y = entropy), color = 'black', size = .3) +
#         geom_line(aes(x = x, y = diff), color = 'darkred', size = .3) +
#         geom_line(
#           data = profile_entropy,
#           aes(x = x, y = entropy),
#           color = 'black',
#           size = .3
#         ) +
#         geom_vline(
#           xintercept = c(lp, rp),
#           color = 'steelblue',
#           linetype = 'dashed'
#         ) +
#         CNAqc:::my_ggplot_theme() +
#         labs(
#           title = bquote("H(x) ~ I = [" * .(round(lp, 2)) * '; ' * .(round(rp, 2)) * ']'),
#           subtitle = paste0(),
#           y = 'Profile',
#           x = bquote('H(x) (black) and' ~ partialdiff * 'H(x) (red)')
#         )
#
#     } # end: with multiple peaks
#
#   # Mutation plots
#   med_coverage = median(snvs_k$DP, na.rm = TRUE)
#
#   mutation_plot = ggplot(snvs_k, aes(VAF, y = ..count.. / sum(..count..))) +
#       geom_histogram(binwidth = 0.01, aes(fill = paste(mutation_multiplicity))) +
#       xlim(0, 1) +
#       CNAqc:::my_ggplot_theme() +
#       scale_color_manual(values = colors) +
#       scale_fill_manual(values = colors) +
#       guides(color = FALSE, fill = guide_legend('Copies')) +
#       labs(
#         y = paste0('Density (simulated ', med_coverage, 'x)'),
#         title = paste0("Mutation multiplicity")
#       )
#
#   if (!(karyotype %in% c('1:1', '1:0')))
#   {
#     d_12 = computation$params$d_12
#     rg_12 = computation$params$rg_12
#     expectation = computation$params$expectation
#
#     mutation_plot = mutation_plot +
#       geom_line(
#         data = data.frame(d_12[1:3]) %>% bind_rows(data.frame(d_12[2 * 1:3])) ,
#         aes(x = x, y = mixture_y),
#         inherit.aes = FALSE,
#         size = .3
#       ) +
#       geom_vline(xintercept = rg_12,
#                  linetype = 3,
#                size = .4) +
#     geom_vline(
#       data = expectation,
#       size = .5,
#       linetype = 'dashed',
#       color = 'indianred2',
#       aes(xintercept = peak, color = label)
#     )
#   }
#
#     # CCF plots
#     CCF_plot = ggplot(snvs_k, aes(CCF, y = ..count.. / sum(..count..))) +
#       geom_histogram(binwidth = 0.01, aes(fill = paste(mutation_multiplicity))) +
#       CNAqc:::my_ggplot_theme() +
#       scale_color_manual(values = colors) +
#       scale_fill_manual(values = colors) +
#       guides(color = FALSE, fill = guide_legend('Copies')) +
#       labs(y = paste0('Density'),
#            title = paste0("CCF for ", karyotype)) +
#       xlim(0, NA)
#
#     figure = cowplot::plot_grid(
#       CCF_plot,
#       mutation_plot,
#       class_plot,
#       entropy_plot,
#       magnitude_plot,
#       # cowplot::plot_grid(magnitude_plot, entropy_plot, align = 'v', nrow = 2),
#       align = 'h',
#       axis = 'b',
#       ncol = 5
#       # scale = c(1, .8)
#     )
#
#     return(figure)
#
# }

mutmult_two_copies_entropy = function(x, karyotype)
{
  cli::cli_rule(
    "Computing mutation multiplicity for karyotype {.field {karyotype}} using the entropy method."
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

  cli::cli_alert_info("Expected Binomial peak(s) for these calls (1 and 2 copies): {.value {expectation$peak}}")

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
  d_1 = CNAqc:::binomial_density(p_1, n, N_bins = 1000)
  d_2 = CNAqc:::binomial_density(p_2, n, N_bins = 1000)

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

  # if(any(duplicated(joint))) joint = joint[!duplicated(joint), ]

  # Prifile the entropy via peak detection
  entropy_profile_x = joint$x
  entropy_profile = joint$entropy

  input_peakdetection = matrix(cbind(x = entropy_profile_x, y = entropy_profile), ncol = 2)
  colnames(input_peakdetection) = c('x', 'y')

  # Peaks detection with these parameters seems to work often
  peaks =  peakPick::peakpick(mat = input_peakdetection,
                              neighlim = 1,
                              deriv.lim = 0.01,
                              peak.min.sd = 0,
                              peak.npos = 1)

  xy_peaks = input_peakdetection[peaks[, 2], , drop = FALSE] %>%
    as_tibble() %>%
    mutate(
      x = x,
      y = x
    )

  xy_peaks

  if(nrow(xy_peaks) == 0) {
    cli::cli_alert_danger("No peaks detected for CCF computation, will not compute values for this karyotype.")

    return(NULL)
  }

  # Points in the centre where there is a violation of the peaks are the actual points we want
  central = entropy_profile_x[which.max(entropy_profile)]

  # signal = entropy_profile
  #
  # J = joint %>%
  #   dplyr::distinct(x, entropy)
  # entropy_profile_x =
  # entropy_profile = joint$entropy
  #
  # dy = J$entropy[-1] - J$entropy[-length(J$entropy)]
  # dx = J$x[-1] - J$x[-length(J$x)]
  #
  # dydx = dy/dx
  #
  #   infl <- c(FALSE, diff(diff(dydx)>0)!=0)
  #   points(J$x[infl ], dydx[infl ], col="blue", pch = 3)
  # abline(v=dydx)
  #
  # mdy = abs(median(dy))
  #
  # lp = rp = which.max(J$entropy)
#
#   repeat{
#     lp = lp - 1
#     if(lp == 1 | abs(dy[lp]) > mdy) break
#   }
#
#   repeat{
#     rp = rp + 1
#     if(rp == length(J$entropy) | abs(dy[rp]) > mdy) break
#   }

  lp = xy_peaks %>% dplyr::filter(x < central) %>% dplyr::arrange(desc(x)) %>% dplyr::filter(row_number() == 1) %>% dplyr::pull(x)
  rp = xy_peaks %>% dplyr::filter(x > central) %>% dplyr::arrange(x) %>% dplyr::filter(row_number() == 1) %>% dplyr::pull(x)

  if(length(lp) == 0 | length(rp) == 0) {
    cli::cli_alert_danger("No suitable range of uncertainty detected for CCF, will not compute values for this karyotype.")

    return(NULL)
  }


  cli::cli_alert_info("Not assignamble area: [{.value {lp}}; {.value {rp}}]")

  # Assignemnts based on lp
  snvs_k = snvs_k %>%
    rowwise() %>%
    mutate(
      mutation_multiplicity = case_when(VAF <= lp ~ '1',
                                        VAF > rp ~ '2',
                                        TRUE ~ 'NA'),
      mutation_multiplicity = as.numeric(mutation_multiplicity),
      CCF = ifelse(
        !is.na(mutation_multiplicity),
        CNAqc:::ccf_adjustment_fun(VAF, B, A, x$purity, mutation_multiplicity),
        NA)
    ) %>%
    ungroup()

  return(list(
    mutations = snvs_k,
    params = list(
      expectation = expectation,
      joint = joint,
      cuts = c(lp, rp),
      method = 'ENTROPY'
    )
  ))
}

mutmult_two_copies_rough = function(x, karyotype)
{
  cli::cli_rule(
    "Computing mutation multiplicity for karyotype {.field {karyotype}} using raw VAF cuts."
  )

  A = as.numeric(strsplit(karyotype, ':')[[1]][1])
  B = as.numeric(strsplit(karyotype, ':')[[1]][2])

  # Karyotype specific mutations - clonal segments
  cl_seg = x$cna %>%
    dplyr::filter(CCF == 1) %>%
    dplyr::pull(segment_id)

  snvs_k = x$snvs %>%
    dplyr::filter(karyotype == !!karyotype, segment_id %in% cl_seg)

  # Expected VAF for 1 and 2 copies of the mutation as for the entropy case
  expectation = CNAqc:::expected_vaf_peak(A, B, x$purity) %>%
    mutate(label = ifelse(mutation_multiplicity == 1, "One copy", "Two copies"))

  med_coverage = median(snvs_k$DP, na.rm = TRUE)

  cli::cli_alert_info("Expected Binomial peak(s) for these calls (1 and 2 copies): {.value {expectation$peak}}.")

  # =-=-=-=-=-=-=-=-=-=-=-
  # Rough-derived heuristic for the detection of points that are difficult to assign
  # =-=-=-=-=-=-=-=-=-=-=-
  p_1 = expectation$peak[1]
  p_2 = expectation$peak[2]

  n = ceiling(med_coverage)

  # We get quantiles as for the entropy
  rg_1 = CNAqc:::binomial_quantile_ranges(p_1, n, quantile_left = 0.01, quantile_right = 0.99)
  rg_2 = CNAqc:::binomial_quantile_ranges(p_2, n, quantile_left = 0.01, quantile_right = 0.99)

  # We create the mixture model: pi_1 * Bin(p1, n) + (1 - pi_1) * Bin(p2, n)
  # as with the entropy
  n_rg_1 = snvs_k %>% filter(VAF > rg_1[1], VAF < rg_1[2]) %>% nrow
  n_rg_2 = snvs_k %>% filter(VAF > rg_2[1], VAF < rg_2[2]) %>% nrow

  # So the algebraic midpoint is (p_2 - p_1)/2, we instead split |p_1-p_2|
  # proportionally to n_rg_1 and n_rg_2, normalised
  mixing = c(n_rg_1, n_rg_2) / (n_rg_1 + n_rg_2)
  t_split = p_1 + (p_2 - p_1) * mixing[1]

  cli::cli_alert_info("Mutations per peak: n = {.value {n_rg_1}}, n = {.value {n_rg_2}}. The hard cut is t = {.value {t_split}}.")

  # Assignemnts based on lp
  snvs_k = snvs_k %>%
    rowwise() %>%
    mutate(
      mutation_multiplicity = case_when(VAF <= t_split ~ '1',
                                        VAF > t_split ~ '2',
                                        TRUE ~ 'NA'),
      mutation_multiplicity = as.numeric(mutation_multiplicity),
      CCF = ifelse(
        !is.na(mutation_multiplicity),
        CNAqc:::ccf_adjustment_fun(VAF, B, A, x$purity, mutation_multiplicity),
        NA)
    ) %>%
    ungroup()


  return(list(
    mutations = snvs_k,
    params = list(
      expectation = expectation,
      cuts = t_split,
      method = 'ROUGH'
    )
  ))
}


plot_mutation_multiplicity_entropy = function(x, karyotype)
{
   # Process
  computation = x$CCF_estimates[[karyotype]]
  snvs_k = computation$mutations

  QC = computation$QC_table$QC

  # Mono-peak
  magnitude_plot = class_plot = entropy_plot = ggplot() + geom_blank()

  colors = RColorBrewer::brewer.pal(n = 3, palette = 'Set1')
  # colors = wesanderson::wes_palette("Royal1", n = 3, type = 'discrete')

  if(!(karyotype %in% c('1:1', '1:0')))
  {
    # Double peaks
    joint = computation$params$joint

    lp = computation$params$cuts[1]
    rp = computation$params$cuts[2]

    # Entropy plot
    entropy_plot = ggplot(joint) +
      geom_line(aes(x = x, y = entropy), color = 'black', size = .3) +
      geom_vline(
        xintercept = c(lp, rp),
        color = 'red',
        linetype = 'dashed',
        size = .3
      ) +
      CNAqc:::my_ggplot_theme() +
      labs(
        title = bquote("NA in [" * .(round(lp, 2)) * '; ' * .(round(rp, 2)) * ']'),
        y = 'Profile',
        x = bquote('Entropy H(x)')
      )

  } # end: with multiple peaks

  # piechart
  ns_counts = sum(is.na(snvs_k$mutation_multiplicity))

  pieplot = snvs_k %>%
    dplyr::group_by(mutation_multiplicity) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::mutate(mutation_multiplicity = paste0(mutation_multiplicity)) %>%
    ggplot(
      aes(x = 1, y = n, fill = mutation_multiplicity)
    ) +
    geom_bar(stat = 'identity') +
    CNAqc:::my_ggplot_theme() +
    scale_fill_manual(values = colors) +
    labs(
      y = paste0(''),
      title = paste0("Assignments (NA = ", ns_counts, ")")
    ) +
    coord_polar(theta = 'y') +
    theme(
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
      ) +
    guides(color = FALSE, fill = guide_legend('Copies'))


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
      y = paste0('Density (', med_coverage, 'x)'),
      title = paste0("Mutation multiplicity")
    )

  if (!(karyotype %in% c('1:1', '1:0')))
  {
    lp = computation$params$cuts[1]
    rp = computation$params$cuts[2]

    # d_12 = computation$params$d_12
    # rg_12 = computation$params$rg_12
    expectation = computation$params$expectation

    mutation_plot = mutation_plot +
      geom_vline(
        data = expectation,
        size = .7,
        linetype = 'dashed',
        color = 'steelblue',
        aes(xintercept = peak, color = label)
      ) +
      geom_vline(
        xintercept = c(lp, rp),
        color = 'red',
        linetype = 'dashed',
        size = .3
      )
    }

  # CCF plots
  CCF_plot = ggplot(snvs_k, aes(CCF, y = ..count.. / sum(..count..))) +
    geom_histogram(bins = 100, aes(fill = paste(mutation_multiplicity))) +
    CNAqc:::my_ggplot_theme() +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    guides(color = FALSE, fill = guide_legend('Copies')) +
    labs(y = paste0('Density'),
         title = paste0("CCF for ", karyotype)) +
    xlim(0, NA)

  # Plot assembly
  # panel = ggpubr::ggarrange(mutation_plot, CCF_plot, pieplot, ncol = 3, nrow = 1, common.legend = T, legend = 'bottom')

  figure = cowplot::plot_grid(
    qc_plot(CCF_plot, QC),
    mutation_plot,
    entropy_plot,
    pieplot,
    align = 'h',
    axis = 'b',
    ncol = 4
  )

  return(figure)
}

plot_mutation_multiplicity_rough = function(x, karyotype)
{
  # Process
  computation = x$CCF_estimates[[karyotype]]
  snvs_k = computation$mutations

  QC = computation$QC_table$QC

  # Mono-peak
  magnitude_plot = class_plot = entropy_plot = ggplot() + geom_blank()

  colors = RColorBrewer::brewer.pal(n = 3, palette = 'Set2')
  # colors = wesanderson::wes_palette("Zissou1", n = 3, type = 'discrete')[c(1,3)]

  # piechart
  ns_counts = sum(is.na(snvs_k$mutation_multiplicity))

  pieplot = snvs_k %>%
    dplyr::group_by(mutation_multiplicity) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::mutate(mutation_multiplicity = paste0(mutation_multiplicity)) %>%
    ggplot(
      aes(x = 1, y = n, fill = mutation_multiplicity)
    ) +
    geom_bar(stat = 'identity') +
    CNAqc:::my_ggplot_theme() +
    scale_fill_manual(values = colors) +
    labs(
      y = paste0(''),
      title = paste0("Assignments (NA = ", ns_counts, ")")
    ) +
    coord_polar(theta = 'y') +
    theme(
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    ) +
    guides(color = FALSE, fill = guide_legend('Copies'))


  # Mutation plots
  med_coverage = median(snvs_k$DP, na.rm = TRUE)
  cuts = computation$params$cuts

  mutation_plot = ggplot(snvs_k, aes(VAF, y = ..count.. / sum(..count..))) +
    geom_histogram(binwidth = 0.01, aes(fill = paste(mutation_multiplicity))) +
    xlim(0, 1) +
    CNAqc:::my_ggplot_theme() +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    guides(color = FALSE, fill = guide_legend('Copies')) +
    labs(
      y = paste0('Density (', med_coverage, 'x)'),
      title = paste0("Mutation multiplicity")
    )
  if (!(karyotype %in% c('1:1', '1:0')))
  {
    lp = computation$params$cuts[1]
    rp = computation$params$cuts[2]

    # d_12 = computation$params$d_12
    # rg_12 = computation$params$rg_12
    expectation = computation$params$expectation

    mutation_plot = mutation_plot +
      geom_vline(
        data = expectation,
        size = .7,
        linetype = 'dashed',
        color = 'steelblue',
        aes(xintercept = peak, color = label)
      ) +
      geom_vline(xintercept = cuts, size = .3, linetype = 'dashed')

  }


  # CCF plots
  CCF_plot = ggplot(snvs_k, aes(CCF, y = ..count.. / sum(..count..))) +
    geom_histogram(bins = 100, aes(fill = paste(mutation_multiplicity))) +
    CNAqc:::my_ggplot_theme() +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    guides(color = FALSE, fill = guide_legend('Copies')) +
    labs(y = paste0('Density'),
         title = paste0("CCF for ", karyotype)) +
    xlim(0, NA)

  # Plot assembly
  # panel = ggpubr::ggarrange(mutation_plot, CCF_plot, pieplot, ncol = 3, nrow = 1, common.legend = T, legend = 'bottom')

  figure = cowplot::plot_grid(
    qc_plot(CCF_plot, QC),
    mutation_plot,
    entropy_plot,
    pieplot,
    align = 'h',
    axis = 'b',
    ncol = 4
  )

  return(figure)
}

plot_mutation_multiplicity_rough_strip = function(x, karyotype)
{
  # Process
  computation = x$CCF_estimates[[karyotype]]
  snvs_k = computation$mutations

  QC = computation$QC_table$QC

  # Mono-peak
  colors = RColorBrewer::brewer.pal(n = 3, palette = 'Set2')
  # colors = wesanderson::wes_palette("Zissou1", n = 3, type = 'discrete')[c(1,3)]

  # Mutation plots
  med_coverage = median(snvs_k$DP, na.rm = TRUE)
  cuts = computation$params$cuts

  mutation_plot = ggplot(snvs_k, aes(VAF, y = ..count.. / sum(..count..))) +
    geom_histogram(binwidth = 0.01, aes(fill = paste(mutation_multiplicity))) +
    xlim(0, 1) +
    CNAqc:::my_ggplot_theme() +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    guides(color = FALSE, fill = guide_legend('Copies')) +
    labs(
      y = paste0('Density (', med_coverage, 'x)'),
      title = paste0("Mutation multiplicity ", karyotype)
    )
  if (!(karyotype %in% c('1:1', '1:0')))
  {
    lp = computation$params$cuts[1]
    rp = computation$params$cuts[2]

    expectation = computation$params$expectation

    mutation_plot = mutation_plot +
      geom_vline(
        data = expectation,
        size = .7,
        linetype = 'dashed',
        color = 'steelblue',
        aes(xintercept = peak, color = label)
      ) +
      geom_vline(xintercept = cuts, size = .3, linetype = 'dashed')

  }

  return(CNAqc:::qc_plot(mutation_plot, QC))
}

plot_mutation_multiplicity_entropy_strip = function(x, karyotype)
{
  # Process
  computation = x$CCF_estimates[[karyotype]]
  snvs_k = computation$mutations

  QC = computation$QC_table$QC

  # Mono-peak
  colors = RColorBrewer::brewer.pal(n = 3, palette = 'Set1')
  # colors = wesanderson::wes_palette("Royal1", n = 3, type = 'discrete')

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
      y = paste0('Density (', med_coverage, 'x)'),
      title = paste0("Mutation multiplicity ", karyotype)
    )

  if (!(karyotype %in% c('1:1', '1:0')))
  {
    lp = computation$params$cuts[1]
    rp = computation$params$cuts[2]

    expectation = computation$params$expectation

    mutation_plot = mutation_plot +
      geom_vline(
        data = expectation,
        size = .7,
        linetype = 'dashed',
        color = 'steelblue',
        aes(xintercept = peak, color = label)
      ) +
      geom_vline(
        xintercept = c(lp, rp),
        color = 'red',
        linetype = 'dashed',
        size = .3
      )
  }

  return(CNAqc:::qc_plot(mutation_plot, QC))
}


# Function to a border to a plot
qc_plot = function(x, QC)
{
  qc = ifelse(QC == "FAIL", "indianred3", 'forestgreen')

  x +
    theme(title = element_text(color = qc),
          panel.border = element_rect(
            colour = qc,
            fill = NA
          ))
}


