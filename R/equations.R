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

  if(purity_from_vaf(vaf,ploidy,multiplicity) > 1){
    
  warning("Incompatible values of VAF, ploidy and multiplicity: the computed score will be unreliable")
    
   }
  
  delta_purity <- (2*multiplicity*delta_vaf)/((multiplicity + vaf*(2-ploidy))^2)
  
  

  return(delta_purity)

}


