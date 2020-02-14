sample_CNA_segments = function(
  karyotypes_p = c(`1:0` = 1, `2:0` = 1, `1:1` = 6, `2:1` = 2, `2:2` = 1),
  breakpoints_rate = 2, 
  ref = "hg19", 
  sex = 'm')
{
  
  K = sum(karyotypes_p)
  karyotypes_p = karyotypes_p/K
  
  if(sum(karyotypes_p) != 1) stop("Sum of karyotypes probabilities must be 1")
  
  reference = CNAqc::chr_coordinates_hg19
  
  if(sex == 'm') reference = reference %>% dplyr::filter(chr != 'chrY')
  if(sex == 'f') reference = reference %>% dplyr::filter(chr != 'chrX')
  
  nchromosomes = nrow(reference)
  
  # Breakpoints are Poisson distributed
  breakpoints_perchr = pio:::nmfy(
    reference$chr,
    rpois(n = nchromosomes, lambda = breakpoints_rate)
  )
  breakpoints_perchr[breakpoints_perchr == 0] = 1
  
  # Chr lengts
  chr_length = pio:::nmfy(
    reference$chr,
    reference$length
  )
  
  calls = lapply(
    reference$chr, 
    function(x)
    {
      nbreaks = breakpoints_perchr[x]
      proportions = MCMCpack::rdirichlet(1, 1:nbreaks) %>% as.vector()
      
      segments_length = round(chr_length[x] * proportions)
      
      karyotypes = sample(
        names(karyotypes_p),
        size = nbreaks,
        prob = karyotypes_p,
        replace = TRUE
      )
      
      Major = strsplit(karyotypes, ':') %>% sapply(FUN = function(x) x[1])
      minor = strsplit(karyotypes, ':') %>% sapply(FUN = function(x) x[2])
      
      df = data.frame(
        chr = x,
        Major = Major %>% as.numeric,
        minor = minor %>% as.numeric,
        from = NA,
        to = NA,
        stringsAsFactors = FALSE
      )
      
      for(i in 1:nrow(df)) {
        df$from[i] = ifelse(i == 1, 0, df$to[i-1])
        df$to[i] = df$from[i] + segments_length[i]
      }
      
      df$length = df$to - df$from
      
      df
    })
  
  Reduce(dplyr::bind_rows, calls) %>%
    dplyr::select(chr, length, from, to, Major, minor) %>%
    tibble::as_tibble()
  
}

sample_VAF_from_CNA = function(
  N,
  segments,
  coverage,
  purity)
{
  # uniformly
  G = segments$length/(sum(segments$length))
  counts = round(G * N)
  
  calls = lapply(
    1:nrow(segments),
    function(i){

      n = counts[i]
      
      if(n == 0) return(data.frame(stringsAsFactors = FALSE))
      
      Major = segments$Major[i]
      minor = segments$minor[i]
      
      copies = 1 + round(runif(n))
      if(Major == 1 & minor == 1) copies = rep(1, n)
      if(Major == 1 & minor == 0) copies = rep(1, n)
      # if(Major == 2 & minor == 0) copies = 1 + round(runif(n))
      # if(Major == 2 & minor == 1) copies = 1 + round(runif(n))
      # if(Major == 2 & minor == 2) copies = 1 + round(runif(n))
      
      vafs = sapply(copies, 
             function(mut.allele)
               CNAqc:::expected_vaf_fun(minor, Major, mut.allele, purity)
      )
      
      lcs = segments$from[i]:segments$to[i]
      lcs = sample(lcs, n)
      
      data.frame(
        chr = segments$chr[i],
        from = lcs,
        to = lcs + 1,
        VAF = vafs,
        ref = sample(c('A', "C", "T", "G"), n, replace = T),
        alt = sample(c('A', "C", "T", "G"), n, replace = T),
        stringsAsFactors = FALSE
      ) %>% tibble::as_tibble()
      
    })
  
  calls = Reduce(dplyr::bind_rows, calls)
  Nc = N
  N = nrow(calls)
  
  
  expected_NV = round(calls$VAF * coverage)
  
  calls$DP = rpois(N, coverage)
  calls$NV = rpois(N, expected_NV)

  # calls$NV = round(calls$VAF * calls$DP)
  calls$VAF = calls$NV/calls$DP
  
  if(nrow(calls) > Nc) calls = calls[1:Nc, ]
  
  calls
}

segments = sample_CNA_segments()
N = 3000
purity = .8
coverage = 120
muts = sample_VAF_from_CNA(N, segments, coverage, purity)

x = CNAqc::init(
  muts, 
  segments,
  purity = purity
)
x

CNAqc::plot_segments(x)
CNAqc::plot_vaf(x)

CNAqc::inspect_segment(x, n = 100)

x = CNAqc::compute_CCF(x)
CNAqc::plot_CCF(x)
