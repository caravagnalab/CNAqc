# Unused
expand_reference_chr_to_arms = function()
{
  p_arm = apply(CNAqc::chr_coordinates_hg19, 1,function(w) {
    data.frame(
      chr = paste0(w['chr'], 'p'),
      chr_original = w['chr'],
      from = w["from"],
      to = w['centromerStart'],
      centromerStart = w['centromerStart'],
      centromerEnd = w['centromerEnd'],
      stringsAsFactors = FALSE
    )
  })

  q_arm = apply(CNAqc::chr_coordinates_hg19, 1,function(w) {
    data.frame(
      chr = paste0(w['chr'], 'q'),
      chr_original = w['chr'],
      from = w["centromerEnd"],
      to = w['to'],
      centromerStart = w['centromerStart'],
      centromerEnd = w['centromerEnd'],
      stringsAsFactors = FALSE
    )
  })

  bind_rows(
    Reduce(bind_rows, p_arm),
    Reduce(bind_rows, q_arm)
  ) %>%
    as_tibble() %>%
    mutate(
      from = as.numeric(from),
      to = as.numeric(to),
      centromerStart = as.numeric(centromerStart),
      centromerEnd = as.numeric(centromerEnd),
      length = to - from
    ) %>%
    select(
      chr, chr_original, length, from, to, centromerStart, centromerEnd
    )

}

split_cna_to_arms = function(clonal_cna)
{
  split_clonal_cna = NULL
  reference_genonme = CNAqc::chr_coordinates_hg19

  for(i in 1:nrow(clonal_cna))
  {
    chr = reference_genonme %>%
      filter(chr == clonal_cna$chr[i])

    # Segment ends before the centromerStart
    if(clonal_cna$to[i] < chr$centromerStart) {
      new = clonal_cna[i, ]
      new$arm = 'p'
      split_clonal_cna = bind_rows(split_clonal_cna, new)
      next
    }

    # Segment start after the centromerEnd
    if(clonal_cna$from[i] >= chr$centromerEnd) {
      new = clonal_cna[i, ]
      new$arm = 'q'
      split_clonal_cna = bind_rows(split_clonal_cna, new)
      next
    }

    # Here it must be split at the centromere, one segment for p, one for q
    new = clonal_cna[i, ]
    new$to = chr$centromerStart
    new$arm = 'p'
    split_clonal_cna = bind_rows(split_clonal_cna, new)

    new = clonal_cna[i, ]
    new$from = chr$centromerEnd
    new$arm = 'q'
    split_clonal_cna = bind_rows(split_clonal_cna, new)
  }

  split_clonal_cna %>%
    mutate(length = to - from) %>%
    select(chr, from, to, length, Major, minor, CCF, arm)
}

# Bayesian test -- this is not used as of today and should be improved to define properly the Bayes Factor
Bayes_test_fragmentation = function(n_short, n_long, integration_range = c(0, .6))
{
  # Bayesian posterior density function for p in a Binomial model with uniform prior
  bayes_posterior_density = function(h, t, r){
    C = factorial(h + t + 1)/(factorial(h) * factorial(t))
    C * {r^h *((1-r)^t)}
  }

  p_domain = seq(0, 1, 0.01)
  density_values = sapply(p_domain, bayes_posterior_density, h = n_short, t = n_long)

  left = integration_range[1]
  right = integration_range[2]

  p_value = integrate(bayes_posterior_density, left, right, h = n_short, t = n_long)
  alternative_p = 1 - right

  plot = ggplot(
    data.frame(y = density_values, x = p_domain)
  ) +
    geom_line(aes(x,y)) +
    geom_vline(xintercept = c(left, right), linetype = 'dashed', color = 'red', size = .3) +
    labs(
      y = "Posterior density",
      x = 'Success probability',
      title = paste0('p = ', round(p_value$value, 6)),
      subtitle = paste0('Segments: S = ', n_short, ", L = ", n_long, '.')
    ) +
    stat_function(fun = bayes_posterior_density,
                  xlim = c(left, right),
                  args = list(h = n_short, t = n_long),
                  geom = "area",
                  fill = 'indianred')


  return(
    list(
      density_values = density_values,
      plot = plot,
      p = p_value$value
    )
  )

}


# Frequentist test
# - two classes, short and long fragments
# - under a uniform belief, if "short" are those shorter than p% of the genome we are testing
#   if assume the Binomial success probability to be just p. Default short is <20% of genome.
frequentist_test_fragmentation = function(n_short,
                                          n_long,
                                          chr,
                                          arm,
                                          testable,
                                          p_cutoff_short = .2,
                                          N_tests = 1,
                                          alpha = 0.01)
{
  if(!testable) return(1)

  # We  consider the probability of observing at least n_short "short" under the null,
  # which constitutes a one-tailed test for whether this "die" is biased towards generating more
  # short segments than expected).
  N = n_long + n_short
  p_value = sum(sapply(n_short:N, function(k) dbinom(k, N, prob = p_cutoff_short)))

  significant = p_value < (alpha/N_tests)

  # Same of this
  # binom.test(n_short, N, p = p_cutoff_short, alternative = 'greater')

  if(significant)
    cli::cli_alert_success(
      "{.field {chr}}{.field {arm}},  p = {.value {p_value}} ~ {.value {N}} segments, {.value {n_short}} short."
      )

  return(p_value)
}

compute_jumps_segments = function(clonal_cna, chr, arm)
{
  clonal_cna = clonal_cna %>%
    filter(chr == !!chr, arm == !!arm) %>%
    mutate(CN_total = Major + minor)

  if(nrow(clonal_cna) <= 1) return(0)

  delta = clonal_cna$CN_total[2:nrow(clonal_cna)] - clonal_cna$CN_total[1:(nrow(clonal_cna) - 1)]
  sum(abs(delta))
}
