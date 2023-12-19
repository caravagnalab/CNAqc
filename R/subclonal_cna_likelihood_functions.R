####################### Clonal Likelihoods #######################################################
### BAF ####
#' Computes expected BAF 
#'
#' @param k1 karyotype (ex '1:1', '1:0'...)
#' @param purity 
#'
#' @return the expected BAF
#' @export
#'
#' @examples
clonal_expected_baf <- function(k1, purity) {
  Ns <- strsplit(k1, ":") %>% unlist() %>% as.numeric()
  nA <- Ns[1]
  nB <- Ns[2]
  
  num <- purity * nB + (1 - purity)
  den <- purity * (nA + nB) + 2 * (1 - purity)
  
  exp_baf <- num / den
  return(exp_baf)
}
#' Computes the Likelihood of the BAF of a set of SNP falling on a single segment
#'
#' @param baf_obs vector with the observed BAF for each bin 
#' @param n number of SNPs / Bin
#' @param k1 karyotype
#' @param purity 
#'
#' @return
#' @export
#'
#' @examples
clonal_baf_ll <- function(baf_obs, n, k1, purity) {
  E_baf <- clonal_expected_baf(k1, purity)
  alpha <- ((n - 2) * E_baf + 1) / (1 - E_baf)
  s <- dbeta(baf_obs, shape1 = alpha, shape2 = n) 
  if (is.infinite(s) && 1000<s){s=1000}
  if (is.infinite(s) && 1000>s){s=-1000}
  return(s)
}

#### DR ####
#' Computes expected depth ratio
#'
#' @param k1 
#' @param purity 
#'
#' @return
#' @export
#'
#' @examples
clonal_expected_dr <- function(k1, purity, ploidy =2){
  na1 <- as.numeric(strsplit(k1, ':')[[1]][1])
  nb1 <- as.numeric(strsplit(k1, ':')[[1]][2])
  
  expected_dpr <- (2 * (1-purity) + purity * (na1 + nb1)) / ploidy
  return(expected_dpr)
}
#' Computes likelihood for the depth ratio of a set of binned SNPs on a single segment 
#'
#' @param dpr_obs vector of DR for each SNP
#' @param n number of SNPs / bin
#' @param k1 karyotype
#' @param purity 
#'
#' @return
#' @export
#'
#' @examples
clonal_dr_ll <- function(dpr_obs, n, k1, purity, ploidy=2){
  expected_dpr <- clonal_expected_dr(k1, purity, ploidy)
  ll <- dgamma(dpr_obs, shape = expected_dpr * sqrt(n) + 1, rate = sqrt(n))
  
  if (is.infinite(ll) && 1000<ll){ll=1000}
  if (is.infinite(ll) && 1000>s){ll=-1000}
  
  return(ll)
}



####################### Subclonal Likelihoods #######################################################
# VAF ####
#' Computes the likelihood for a single SNV
#'
#' @param k number of reads in the region carrying the variant (NV)
#' @param n number of reads in the region (DP)
#' @param ps vector of the VAF of the expected peaks (ex. c(.5) for a clonal diploid region, with purity 1)
#'
#' @return
#' @export
#'
#' @examples
vaf_ll <- function(k, n, ps){
  #AVERAGE of the LIKELIHOOD
  s <- 0
  for (peak in 1:length(ps)){
    #s <- s + max(dbinom(k, size=n, prob=ps[peak]))
    s <- s + (1/length(ps)) * dbinom(k, size=n, prob=ps[peak])
    if (is.infinite(s) && 1000<s){s=1000}
    if (is.infinite(s) && 1000>s){s=-1000}
  }
  return(s)
}
#' Computes the log-likelihood for a set of SNVs falling on a single segment 
#'
#' @param NV vector of number of reads in the region carrying the variants (NV)
#' @param DP vector of number of reads in the region (DP)
#' @param peaks vector of the VAF of the expected peaks (ex. c(.5) for a clonal diploid region, with purity 1)
#'
#' @return the log-likelihood of the set of SNVs
#' @export
#'
#' @examples
VAF_LL <- function(NV, DP, peaks){
  vaf_ll <- 0
  
  #peaks = peaks[peaks>.15]
  
  for (i in 1:length(NV)) {
    v <- vaf_ll(k = NV[i], n = DP[i], ps = peaks$peak) #peaks$peak
    vaf_ll <- vaf_ll + log(v)
  }
  if (is.infinite(vaf_ll) && 100<vaf_ll){vaf_ll=1000000000}
  if (is.infinite(vaf_ll) && 100>vaf_ll){vaf_ll=-100000000}
  return(vaf_ll)
}

# BAF ####
#' Computes expected BAF in the case of two subclones
#'
#' @param k1 karyorype subclone 1
#' @param k2 karyorype subclone 2
#' @param purity 
#' @param ccf 
#' @param g1 genotype subclone 1 (the mathernal and paternal alleles are named A and B, plus a number indicating when they were originated : A1B1 = 1:1 )
#' @param g2 genotype subclone 2
#'
#' @return expected BAF
#' @export
#'
#' @examples
expected_baf <- function(k1, k2, purity, ccf, g1, g2){
  nA1 <- lengths(regmatches(g1, gregexpr("A", g1)))
  nB1 <- lengths(regmatches(g1, gregexpr("B", g1)))
  nA2 <- lengths(regmatches(g2, gregexpr("A", g2)))
  nB2 <- lengths(regmatches(g2, gregexpr("B", g2)))
  
  
  num <- min(c( (nA1 * ccf + nA2 * (1 - ccf))*purity , (nB1 * ccf + nB2 * (1 - ccf))*purity)) + (1- purity)
  den <- purity * ((nA1 + nB1) * ccf + (1 - ccf) * (nA2 + nB2)) + 2 * (1 - purity)
  exp_baf <- num / den
  return(exp_baf)
}
#' Computes the likelihood of the BAF of a set of binned SNP sitting on a subclonal segment
#'
#' @param baf_obs vector of observed BAF for each bin
#' @param n number SNPs/bin
#' @param k1 karytype subclone 1
#' @param k2 karytype subclone 2
#' @param purity 
#' @param ccf 
#' @param g1 genotype subclone 1
#' @param g2 genotype subclone 2
#'
#' @return likelihood of the BAF of a set of binned SNP sitting on a subclonal segment
#' @export
#'
#' @examples
BAF_LL <- function(baf_obs, n, k1, k2, purity, ccf, g1, g2){
  E_baf <- expected_baf(k1, k2, purity, ccf, g1, g2)
  alpha <- ((n - 2) * E_baf + 1) / (1 - E_baf)
  s <- dbeta(baf_obs, shape1 = alpha, shape2 = n) 
  if (is.infinite(s) && 1000<s){s=1000}
  if (is.infinite(s) && 1000>s){s=-1000}
  return(s)
}


# DR ####
#' Computes expected DR in the case of two subclones
#'
#' @param k1 karytype subclone 1
#' @param k2 karytype subclone 2
#' @param purity 
#' @param ccf 
#' @param ploidy average ploidy of the sample, default = 2
#'
#' @return expected DR in the case of two subclones
#' @export
#'
#' @examples
expected_dr <- function(k1, k2, purity, ccf, ploidy = 2){
  na1 <- as.numeric(strsplit(k1, ':')[[1]][1])
  nb1 <- as.numeric(strsplit(k1, ':')[[1]][2])
  na2 <- as.numeric(strsplit(k2, ':')[[1]][1])
  nb2 <- as.numeric(strsplit(k2, ':')[[1]][2])
  expected_dpr <- (2*(1-purity) + purity*(ccf*(na1 + nb1) + (1- ccf)*(na2 + nb2))) / ploidy 
  return(expected_dpr)
  
}
#' Computes the likelihood of the DR of a set of binned SNP sitting on a subclonal segment
#'
#' @param dpr_obs vector of observed DR for each bin
#' @param n number of SNPs/bin
#' @param k1 karytype subclone 1
#' @param k2 karytype subclone 2
#' @param purity 
#' @param ccf 
#' @param ploidy average ploidy of the sample, default = 2
#'
#' @return likelihood of the DR of a set of binned SNP sitting on a subclonal segment
#' @export
#'
#' @examples
DR_LL <- function(dpr_obs, n, k1, k2, purity, ccf, ploidy = 2){
  expected_dpr <- expected_dr(k1, k2, purity, ccf, ploidy)
  ll <- dgamma(dpr_obs, shape = expected_dpr * sqrt(n) + 1, rate = sqrt(n))
  if (is.infinite(ll) && 1000<ll){ll=1000}
  if (is.infinite(ll) && 1000>ll){ll=-1000}
  return(ll)
}







