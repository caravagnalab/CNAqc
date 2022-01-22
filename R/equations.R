# Expected VAF for a mutation mapping to a clonal CNA segment
# with certain minor/ Major alelles, with a number of copies (of the mutation)
# fixed (either 1 or more), with a sample purity.
# m - minor allele
# M - Major allele
# p - purity
# mut.allele - mutation multiplicity
expected_vaf_fun = function(m, M, mut.allele, p)
{
  P = m + M

  expected_mutant_reads = mut.allele * p
  expected_sequencing_depth = 2 * (1 - p) + p * P
  expected_mutant_reads / expected_sequencing_depth
}

# Expected VAF for a general peak, using expected_vaf_fun
expectations_generalised = function(m, M, p, karyotype = NULL)
{
  if(!is.null(karyotype))
  {
    karyotype = strsplit(karyotype, ':')[[1]]
    m = karyotype[2] %>% as.numeric()
    M = karyotype[1] %>% as.numeric()
  }

  minor_peaks = lapply(1:m, function(i){
    data.frame(
      minor = m,
      Major = M,
      ploidy = m + M,
      multiplicity = i,
      purity = p,
      peak = expected_vaf_fun(m, M, mut.allele = i, p)
    )
  })

  Major_peaks = lapply(1:M, function(i){
    data.frame(
      minor = m,
      Major = M,
      ploidy = m + M,
      multiplicity = i,
      purity = p,
      peak = expected_vaf_fun(m, M, mut.allele = i, p)
    )
  })

  return(
    bind_rows(
      Reduce(bind_rows, minor_peaks),
      Reduce(bind_rows, Major_peaks)
    ) %>%
      distinct() %>%
      mutate(karyotype = paste(Major, minor, sep = ':')) %>%
      filter(multiplicity > 0)
  )
}

# Expectations for subclonal peaks - linear evolution model
expectations_subclonal_linear = function(CCF_1, karyotype_1, karyotype_2, purity)
{
  CCF_2 = 1 - CCF_1
  ploidy_1 = strsplit(karyotype_1, ":")[[1]] %>% as.numeric() %>% sum
  ploidy_2 = strsplit(karyotype_2, ":")[[1]] %>% as.numeric() %>% sum

  peak1 = peak2 = peak3 = peak4 = peak5 = NA
  role1 = role2 = role3 = role4 = role5 = NA

  # Equations - denominator is always the same
  den = 2 * (1 - purity) + ploidy_1 * CCF_1 * purity + ploidy_2 * CCF_2 *
    purity

  if ((karyotype_1 == "1:0" & karyotype_2 == "1:1") |
      (karyotype_1 == "1:1" & karyotype_2 == "1:0"))
  {
    # 1:0 1:1 - 1:1 1:0
    peak1 = ((1 * CCF_1 + 1 * CCF_2) * purity) / den
    peak2 = ((1 * CCF_1) * purity) / den
    peak3 = ((1 * CCF_2) * purity) / den

    role1 = 'shared'
    role2 = role3 = 'private'
  }

  if (karyotype_1 == "1:0" & karyotype_2 == "2:1")
  {
    # 1:0 2:1
    peak1 = ((1 * CCF_1 + 2 * CCF_2) * purity) / den
    peak2 = ((1 * CCF_1 + 1 * CCF_2) * purity) / den
    peak3 = ((1 * CCF_1) * purity) / den
    peak4 = ((1 * CCF_2) * purity) / den

    role1 = role2 = 'shared'
    role3 = role4 = 'private'
  }

  if (karyotype_1 == "2:1" & karyotype_2 == "1:0")
  {
    # 2:1 1:0
    peak1 = ((2 * CCF_1 + 1 * CCF_2) * purity) / den
    peak2 = ((1 * CCF_1 + 1 * CCF_2) * purity) / den
    peak3 = ((1 * CCF_1) * purity) / den
    peak4 = ((1 * CCF_2) * purity) / den

    role1 = role2 = 'shared'
    role3 = role4 = 'private'
  }

  
  if (karyotype_1 == "1:0" & karyotype_2 == "2:2")
  {
     # 1:0 2:2 
    peak1 = ((1 * CCF_1 + 2 * CCF_2) * purity) / den
    peak2 = ((1 * CCF_1) * purity) / den
    peak3 = ((1 * CCF_2) * purity) / den
    peak4 = ((2 * CCF_2) * purity) / den

    role1 = 'shared'
    role2 = role3 = role4 = 'private'
  }
  
  if (karyotype_1 == "2:2" & karyotype_2 == "1:0")
  {
     # 2:2 1:0
    peak1 = ((2 * CCF_1 + 1 * CCF_2) * purity) / den
    peak2 = ((1 * CCF_1) * purity) / den
    peak3 = ((1 * CCF_2) * purity) / den
    peak4 = ((2 * CCF_1) * purity) / den

    role1 = 'shared'
    role2 = role3 = role4 = 'private'
  }


  if (karyotype_1 == "1:1" & karyotype_2 == "2:1")
  {
    # 1:1 2:1
    peak1 = ((1 * CCF_1 + 2 * CCF_2) * purity) / den
    peak2 = ((1 * CCF_1) * purity) / den
    peak3 = ((1 * CCF_2) * purity) / den

    role1 = 'shared'
    role2 = role3 = 'private'
  }
  
  if (karyotype_1 == "1:0" & karyotype_2 == "2:0")
  {
    # 1:0 2:0
    peak1 = ((1 * CCF_1 + 2 * CCF_2) * purity) / den
    peak2 = ((1 * CCF_1) * purity) / den
    peak3 = ((2 * CCF_1) * purity) / den
    peak4 = ((1 * CCF_2) * purity) / den

    role1 = 'shared'
    role2 = role3 = role4 = 'private'
  }
  
  if (karyotype_1 == "2:0" & karyotype_2 == "1:0")
  {
    # 2:0 1:0
    peak1 = ((2 * CCF_1 + 1 * CCF_2) * purity) / den
    peak2 = ((1 * CCF_1) * purity) / den
    peak3 = ((2 * CCF_2) * purity) / den
    peak4 = ((1 * CCF_2) * purity) / den

    role1 = 'shared'
    role2 = role3 = role4 = 'private'
  }

  if (karyotype_1 == "2:1" & karyotype_2 == "1:1")
  {
    # 2:1 1:1
    peak1 = ((1 * CCF_2 + 2 * CCF_1) * purity) / den
    peak2 = ((1 * CCF_1) * purity) / den
    peak3 = ((1 * CCF_2) * purity) / den

    role1 = 'shared'
    role2 = role3 = 'private'
  }

  if (karyotype_1 == "1:1" & karyotype_2 == "2:0")
  {
    # 1:1 2:0
    peak1 = ((1 * CCF_1 + 2 * CCF_2) * purity) / den
    peak2 = ((1 * CCF_1) * purity) / den
    peak3 = ((1 * CCF_2) * purity) / den

    role1 = 'shared'
    role2 = role3 = 'private'
  }

  if (karyotype_1 == "2:0" & karyotype_2 == "1:1")
  {
    # 2:0 1:1
    peak1 = ((1 * CCF_2 + 2 * CCF_1) * purity) / den
    peak2 = ((1 * CCF_1) * purity) / den
    peak3 = ((1 * CCF_2) * purity) / den

    role1 = 'shared'
    role2 = role3 = 'private'
  }
  
  if (karyotype_1 == "1:1" & karyotype_2 == "2:2")
  {
    # 1:1 2:2
    peak1 = ((1 * CCF_1 + 2 * CCF_2) * purity) / den
    peak2 = ((1 * CCF_1) * purity) / den
    peak3 = ((1 * CCF_2) * purity) / den
    peak3 = ((2 * CCF_2) * purity) / den

    role1 = 'shared'
    role2 = role3 = role4= 'private'
  }
  
   if (karyotype_1 == "2:2" & karyotype_2 == "1:1")
  {
    # 2:2 1:1
    peak1 = ((2 * CCF_1 + 1 * CCF_2) * purity) / den
    peak2 = ((1 * CCF_1) * purity) / den
    peak3 = ((1 * CCF_2) * purity) / den
    peak3 = ((2 * CCF_1) * purity) / den

    role1 = 'shared'
    role2 = role3 = role4= 'private'
  }

  if (karyotype_1 == "2:1" & karyotype_2 == "2:2")
  {
    # 2:1 2:2
    peak1 = ((2 * CCF_1 + 2 * CCF_2) * purity) / den
    peak2 = ((1 * CCF_1 + 1 * CCF_2) * purity) / den
    peak3 = ((1 * CCF_1 + 2 * CCF_2) * purity) / den
    peak4 = ((1 * CCF_1) * purity) / den
    peak5 = ((1 * CCF_2) * purity) / den

    role1 = role2 = role3 = 'shared'
    role4 = role5 = 'private'
  }

  if (karyotype_1 == "2:2" & karyotype_2 == "2:1")
  {
    # 2:2 2:1
    peak1 = ((2 * CCF_1 + 2 * CCF_2) * purity) / den
    peak2 = ((1 * CCF_1 + 1 * CCF_2) * purity) / den
    peak3 = ((1 * CCF_2 + 2 * CCF_1) * purity) / den
    peak4 = ((1 * CCF_1) * purity) / den
    peak5 = ((1 * CCF_2) * purity) / den

    role1 = role2 = role3 = 'shared'
    role4 = role5 = 'private'
  }

  if ((karyotype_1 == "2:1" &
       karyotype_2 == "2:0") |
      (karyotype_1 == "2:0" & karyotype_2 == "2:1"))
  {
    # 2:1 2:0 - 2:0 2:1
    peak1 = ((2 * CCF_1 + 2 * CCF_2) * purity) / den
    peak2 = ((1 * CCF_1 + 1 * CCF_2) * purity) / den
    peak3 = ((1 * CCF_1) * purity) / den
    peak4 = ((1 * CCF_2) * purity) / den

    role1 = role2 = 'shared'
    role3 = role4 = 'private'
  }

  if ((karyotype_1 == "2:2" &
       karyotype_2 == "2:0") |
      (karyotype_1 == "2:0" & karyotype_2 == "2:2"))
  {
    # 2:2, 2:0 - 2:0 2:2
    peak1 = ((2 * CCF_1 + 2 * CCF_2) * purity) / den
    peak2 = ((1 * CCF_1 + 1 * CCF_2) * purity) / den
    peak3 = ((1 * CCF_1) * purity) / den
    peak4 = ((1 * CCF_2) * purity) / den

    role1 = role2 = 'shared'
    role3 = role4 = 'private'
  }

  df = data.frame(
    karyotype_1 = karyotype_1,
    karyotype_2 = karyotype_2,
    CCF_1 = CCF_1,
    CCF_2 = CCF_2,
    peak = c(peak1, peak2, peak3, peak4, peak5),
    role = c(role1, role2, role3, role4, role5)
  )

  df = df[complete.cases(df), ]

  return(df)
}

# Expectations for subclonal peaks - branching evolution model
expectations_subclonal_branching = function(CCF_1, karyotype_1, karyotype_2, purity)
{
  CCF_2 = 1 - CCF_1

  # Subclones ploidies
  ploidy_1 = strsplit(karyotype_1, ':')[[1]] %>% as.numeric() %>% sum # p_1
  ploidy_2 = strsplit(karyotype_2, ':')[[1]] %>% as.numeric() %>% sum # p_2

  # Multiplicities expected - m_i_j ~ multiplicity j in clone i
  m_1_1 = m_1_2 = 1
  m_2_1 = m_2_2 =  2

  # Common denominator
  denominator =
    (1 - purity) * 2 +
    ploidy_1 * CCF_1 * purity +
    ploidy_2 * CCF_2 * purity

  # Numerator subclones -  v_i_j ~ reads from peak at multiplicity j in clone i
  v_1_1 = m_1_1 * CCF_1 * purity
  v_2_1 = m_2_1 * CCF_1 * purity

  v_1_2 = m_1_2 * CCF_2 * purity
  v_2_2 = m_2_2 * CCF_2 * purity

  # VAF peaks subclones ~ peak at multiplicity j in clone i
  peak_1_1 = v_1_1 / denominator
  peak_2_1 = v_2_1 / denominator

  peak_1_2 = v_1_2 / denominator
  peak_2_2 = v_2_2 / denominator

  # Filter unused peaks
  if (karyotype_1 %in% c("1:0", "1:1"))
    peak_2_1 = NA
  if (karyotype_2 %in% c("1:0", "1:1"))
    peak_2_2 = NA

  df = data.frame(
    karyotype_1 = karyotype_1,
    karyotype_2 = karyotype_2,
    CCF_1 = CCF_1,
    CCF_2 = CCF_2,
    peak = c(peak_1_1, peak_2_1, peak_1_2, peak_2_2),
    role = 'private'
  )

  df[complete.cases(df), ]
}

# Expectations for subclonal peaks - all models
expectations_subclonal = function(CCF_1, karyotype_1, karyotype_2, purity)
{
  el = expectations_subclonal_linear(CCF_1, karyotype_1, karyotype_2, purity)
  if(nrow(el) > 0) el = el %>% mutate(model = 'linear')

  es = expectations_subclonal_branching(CCF_1, karyotype_1, karyotype_2, purity)
  if(nrow(es) > 0) es = es %>% mutate(model = 'branching')

  bind_rows(el, es)
}

# Compute CCF values from mutation multiplicity
# m - minor allele
# M - Major allele
# p - purity
# mut.allele - mutation multiplicity
ccf_adjustment_fun = function(v, m, M, p, mut.allele = 1)
{
  CN = as.numeric(m) + as.numeric(M)
  v = as.numeric(v)
  p = as.numeric(p)
  mut.allele = as.numeric(mut.allele)

  v * ((CN - 2) * p + 2) / (mut.allele * p)
}


# Compute sample purity from mutation multiplicity
# m - minor allele
# M - Major allele
# p - purity
# mut.allele - mutation multiplicity
purity_estimation_fun = function(v, m, M, mut.allele = 1)
{
  CN = as.numeric(m) + as.numeric(M)
  v = as.numeric(v)
  mut.allele = as.numeric(mut.allele)

  (2 * v) / (mut.allele + v * (2 - CN))
}

# Tetraploid (m = M = 2), VAF 50% -- pure tumour!
# purity_estimation_fun(v = .5, m = 2, M = 2, mut.allele = 2)
#
# # Triploid (m = 1, M = 2), VAF 2/3%  -- pure tumour!
# purity_estimation_fun(v = 2/3, m = 1, M = 2, mut.allele = 2)
#
# purity_estimation_fun(v = .66, m = 1, M = 2, mut.allele = 2)
#
# # Diploid balanced (m = M = 1), VAF 50% -- pure tumour!
# purity_estimation_fun(v = .5, m = 1, M = 1, mut.allele = 1)
# purity_estimation_fun(v = .06, m = 1, M = 1, mut.allele = 1)
# purity_estimation_fun(v = .06, m = 2, M = 2, mut.allele = 2)
#
# purity_estimation_fun(v = .06, m = 0, M = 1, mut.allele = 1)
#
# purity_estimation_fun(v = .06, m = 1, M = 1, mut.allele = 1)

# Compute VAF values from CCF and mutation multiplicity
# m - minor allele
# M - Major allele
# p - purity
# mut.allele - mutation multiplicity
vaf_from_ccf = function(ccf, m, M, p, mut.allele = 1)
{
  CN = as.numeric(m) + as.numeric(M)
  ccf = as.numeric(ccf)
  p = as.numeric(p)
  mut.allele = as.numeric(mut.allele)

  (mut.allele * p * ccf) / ((CN - 2) * p + 2)
}

# vaf_from_ccf(1, 1, 1, 1, 1)
# vaf_from_ccf(.3, 1, 1, 1, 1)
# vaf_from_ccf(1, 1, 2, 1, 1)

# formula: delta_vaf= 2*multiplicty*epsilon_error/((2+purity(ploidy-2))^2)
#compute delta_vaf for all karyotypes and multiplicity given epsilon_error and purity
delta_vaf_karyo = function(epsilon_error, purity)
{
  compute_vaf_error <-
    function(epsilon_error,
             purity,
             multiplicity,
             ploidy) {
      delta_vaf = (2 * multiplicity * epsilon_error) / ((2 + purity * (ploidy -
                                                                         2)) ^ 2)
      return(delta_vaf)
    }

  delta_10 <-
    tibble(
      karyotype = "1:0",
      multiplicity = 1,
      delta_vaf = compute_vaf_error(epsilon_error, purity, 1, 1)
    )

  delta_11 <-
    tibble(
      karyotype = "1:1",
      multiplicity = 1,
      delta_vaf = compute_vaf_error(epsilon_error, purity, 1, 2)
    )

  delta_20 <- tibble(
    karyotype = "2:0",
    multiplicity = c(1, 2),
    delta_vaf = c(
      compute_vaf_error(epsilon_error, purity, 1, 2),
      compute_vaf_error(epsilon_error, purity, 2, 2)
    )
  )

  delta_21 <- tibble(
    karyotype = "2:1",
    multiplicity = c(1, 2),
    delta_vaf = c(
      compute_vaf_error(epsilon_error, purity, 1, 3),
      compute_vaf_error(epsilon_error, purity, 2, 3)
    )
  )

  delta_22 <- tibble(
    karyotype = "2:2",
    multiplicity = c(1, 2),
    delta_vaf = c(
      compute_vaf_error(epsilon_error, purity, 1, 4),
      compute_vaf_error(epsilon_error, purity, 2, 4)
    )
  )

  Delta_vaf <- rbind(delta_10, delta_11, delta_20, delta_21, delta_22)

  return(Delta_vaf)
}


# invert vaf(purity) and compute purity from vaf,ploidy and muliplicity
purity_from_vaf <- function(vaf,ploidy,multiplicity){

  purity <- (2*vaf)/(multiplicity + (2-ploidy)*vaf)

  return(purity)

}


# get vaf,delta_vaf and karyotype(ploidy,multiplicty) and return delta_purity
compute_delta_purity <- function(vaf,delta_vaf,ploidy,multiplicity){

  if (purity_from_vaf(vaf, ploidy, multiplicity) > 1) {
    warning(
      "Incompatible VAF, ploidy and multiplicity: the computed score might be unreliable"
    )
  }

  delta_purity <- (2*multiplicity*delta_vaf)/((multiplicity + vaf*(2-ploidy))^2)



  return(delta_purity)

}


