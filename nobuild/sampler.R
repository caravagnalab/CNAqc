#' Sample an allele-specific CNA profile.
#'
#' @description Returns a chromosome-level segmentation of a genome
#' (hg19 or c) with associated allele-specific clonal Copy Number
#' Alterations. The karyotypes probabilities (multinomials) are given in
#' input, together with a Poisson rate parameter to generate breakpoints.
#'
#' The function returns a tibble reporting information about the
#' sampled copy number segments.
#'
#' @param karyotypes_p A named vector of multinomial coefficients
#' in the form \code{`M:m` = p} denoting the un-normalised probability
#' mass of generating a segment with \code{"M"} copies of the major allele,
#' and \code{"m"} copies of the minor allele.
#' @param breakpoints_rate Poisson rate $\lambda>0$ to sample the number of
#' breakpoints per chromosome. The sampling process ensures that
#'
#' @param reference_genome
#' @param sex
#'
#' @return
#' @export
#'
#' @examples
sample_CNA_segments = function(
  karyotypes_p = c(`1:0` = 1, `2:0` = 1, `1:1` = 6, `2:1` = 2, `2:2` = 1),
  breakpoints_rate = 2,
  reference_genome = "GRCh38",
  sex = 'm')
{
  if(length(karyotypes_p) == 0) stop("No input karyotypes.")
  if(any(karyotypes_p < 0)) stop("The probability per karyotype must be positive.")
  if(breakpoints_rate < 0) stop("Breakpoints rate must be positive.")

  cli::cli_h1("Allele-specific CNAs samples")

  cli::cli_h2("Karyotypes ~ Multinomial({.field {karyotypes_p}})")

  # Normalise karyotypes probability
  K = sum(karyotypes_p)
  karyotypes_p = karyotypes_p/K

  for(k in karyotypes_p %>% sort(decreasing = TRUE) %>% names)
    cli::cli_alert("Karyotype {.field {k}} with p = {.value {karyotypes_p[k]}}")

  # Get reference
  reference = NULL
  if(reference_genome == 'hg19' | reference_genome == "GRCh37")
    reference = CNAqc::chr_coordinates_hg19

  if(reference_genome == "GRCh38")
    reference = CNAqc::chr_coordinates_GRCh38

  # if(sex == 'm') reference = reference %>% dplyr::filter(chr != 'chrY')
  if(sex == 'f') reference = reference %>% dplyr::filter(chr != 'chrY')

  nchromosomes = nrow(reference)

  # Breakpoints are Poisson distributed
  breakpoints_perchr = pio:::nmfy(
    reference$chr,
    rpois(n = nchromosomes, lambda = breakpoints_rate)
  )
  # breakpoints_perchr[breakpoints_perchr == 0] = 1

  cli::cli_h2("Breakpoints per chromosome: B ~ Poisson({.field {breakpoints_rate}})")

  # breakpoints_perchr %>%
  #   as_tibble() %>%
  #   mutate(chr = breakpoints_perchr %>% names) %>%
  #   rename(breakpoints = value) %>%
  #   arrange(breakpoints %>% desc) %>%
  #   print

  for(br in breakpoints_perchr %>% names)
    cat(
      sprintf("%5s", br),
      rep("\u25A0", breakpoints_perchr[br]) %>% paste(collapse = ''),
      '\n')

  chr_length = pio:::nmfy(
    reference$chr,
    reference$length
  )

  # Sample breakpoints per chromosome
  calls = lapply(
    reference$chr,
    function(x)
    {
      # Sample the length of the segments
      nbreaks = breakpoints_perchr[x]

      # Special case, no breakpoints
      if(nbreaks == 0)
      {
        sgm = sample(
          names(karyotypes_p),
          size = 1,
          prob = karyotypes_p
        )

        df = data.frame(
          chr = x,
          Major = strsplit(sgm, ':') %>% sapply(FUN = function(x) x[1]) %>% as.numeric,
          minor = strsplit(sgm, ':') %>% sapply(FUN = function(x) x[2]) %>% as.numeric,
          from = 0,
          to = chr_length[x],
          stringsAsFactors = FALSE
        )
        df$length = df$to - df$from

        return(df)
      }

      # General case

      proportions = runif(nbreaks) %>% as.vector()
      proportions = proportions/sum(proportions)

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

  cat("\n")
  cli::cli_h2("Segments")

  cnas = Reduce(dplyr::bind_rows, calls) %>%
    dplyr::select(chr, from, to, length, Major, minor) %>%
    tibble::as_tibble() %>%
    mutate(ploidy = Major + minor)

  print(cnas)

  cnas %>% return
}

#' Sample clonal SNVs matching CNA segments.
#'
#' @description Given an allele-specific segmentation of the
#' tumour genome, this function returns a tibble of mutations
#' associated to the input clonal copy number segments.
#'
#' First, it assigns somatic mutations
#' to each segment proportionally to the ploidy, then determines
#' the multiplicitly of each mutation (randomly) and computes
#' the associated theoretical VAF.
#'
#' Finally, it samples Poisson-distributed coverage and then
#' simulate read counts data per mutation.
#'
#' @param N Number of required mutations.
#' @param segments Input clonal allele-specific segmentation.
#' @param coverage Mean coverage per allele.
#' @param purity Sample purity.
#'
#' @return
#' @export
#'
#' @examples
sample_clonal_mutations_from_CNA = function(
  N,
  segments,
  coverage,
  purity)
{
  cli::cli_h1("Clonal mutation sampling from allele-specific CNAs")

  cli::cli_h2("Mapping mutations consistently with genome ploidy")

  # proportionally to size and ploidy
  G = segments$length * segments$ploidy
  G = G / sum(G)
  counts = round(G * N)

  cli::cli_h2("Generating read-counts data, per-allele coverage ~ Poisson({.field {coverage}})")

  calls = lapply(
    1:nrow(segments),
    function(i){

      n = counts[i]
      ploidy = segments$ploidy[i]

      # Adjust coverage as needed
      s_coverage = coverage * ploidy

      if(n == 0) return(data.frame(stringsAsFactors = FALSE))

      Major = segments$Major[i]
      minor = segments$minor[i]

      copies = 1 + round(runif(n))
      if(Major == 1 & minor == 1) copies = rep(1, n)
      if(Major == 1 & minor == 0) copies = rep(1, n)

      # True VAF
      vafs = sapply(copies,
             function(mut.allele)
               CNAqc:::expected_vaf_fun(minor, Major, mut.allele, purity)
      )

      lcs = segments$from[i]:segments$to[i]
      lcs = sample(lcs, n)

      df = data.frame(
        chr = segments$chr[i],
        from = lcs,
        to = lcs + 1,
        ref = sample(c('A', "C", "T", "G"), n, replace = T),
        alt = sample(c('A', "C", "T", "G"), n, replace = T),
        multiplicity = copies,
        Major = segments$Major[i],
        minor = segments$minor[i],
        ploidy = segments$ploidy[i],
        th_VAF = vafs,
        stringsAsFactors = FALSE
      ) %>%
        as_tibble() %>%
        arrange(from)

      df$DP = rpois(nrow(df), s_coverage)
      df$NV = rbinom(nrow(df), size = df$DP, prob = df$th_VAF)
      df$VAF = df$NV/df$DP

      return(df)
    })

  calls = Reduce(dplyr::bind_rows, calls)

  print(calls)

  calls
}


###############
# Example run #
###############

# Input segmentation
segments = sample_CNA_segments(
  breakpoints_rate = 2
  )

# 3000 clonal mutations, 80% purity, per-allele coverage 50 (means
# 100x coverage on a diploid area)
N = 3000
purity = .8
coverage = 50
muts = sample_clonal_mutations_from_CNA(N, segments, coverage, purity)

# Test the sampled data
x = CNAqc::init(
  muts,
  segments,
  purity = purity
)

CNAqc::plot_segments(x)
CNAqc::plot_data_histogram(x)
CNAqc::plot_data_histogram(x, which = 'DP')

# Check peaks
x = CNAqc::analyze_peaks(x, min_karyotype_size = .05)

CNAqc::plot_peaks_analysis(x, empty_plot = FALSE)

CNAqc::inspect_segment(x, n = 100)

# CCF
x = CNAqc::compute_CCF(x)
CNAqc::plot_CCF(x, empty_plot = FALSE)


