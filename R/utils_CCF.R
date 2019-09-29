# Density values from a Binmoial distribution with number of trials `n`
# and success probability `p`. Returns a data.frame with the density
binomial_density = function(p, n)
{
  # Domain
  domain = seq(0, 1, 0.01)

  # Bernoulli succesfull trials
  NV = round(domain * n)

  data.frame(
    x = NV/n,
    y = dbinom(NV, size = n, prob = p)
  )
}

# Quantile from a Binmoial distribution with number of trials `n`
# and success probability `p`. Returns a data.frame with the density
binomial_quantile_ranges = function(p, n, quantile_left = 0.01, quantile_right = 1 - quantile_left)
{
  c(
    qbinom(quantile_left, n, p)/n,
    qbinom(quantile_right, n, p)/n
  )
}

# Create a entropy profile for a 2-class model
entropy_profile_2_class = function(d_1, d_2)
{
  joint = data.frame(
    A = d_1$mixture_y,
    B = d_2$mixture_y
  )
  C = apply(joint, 1, sum)

  joint$A = joint$A/C
  joint$B = joint$B/C

  # Entropy of the mixture
  joint$x = d_1$x
  joint$entropy = - apply(joint * log(joint), 1, sum)
  joint$entropy[is.na(joint$entropy)] = 0

  # Absolute derivative (magnitude) of the entropy
  sz = length(joint$entropy)
  joint$diff = abs(c(joint$entropy[2:sz] - joint$entropy[1:(sz - 1)], 0))

  joint
}
