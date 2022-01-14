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

# data('example_dataset_CNAqc', package = 'CNAqc')
# x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna,example_dataset_CNAqc$purity)

analyze_peaks_general = function(x, n_min = 50, epsilon = 0.03,
                                 kernel_adjust = 1, n_bootstrap = 1)
{
  aux_peak_calling_general = function(mutations){
    xy_peaks = den = NULL

    # Smoothed Gaussian kernel for VAF
    y = mutations %>% dplyr::pull(VAF)

    den = density(y, kernel = 'gaussian', adjust = kernel_adjust, na.rm = T)
    in_range = den$x >= min(y, na.rm = T) & den$x <= max(y, na.rm = T)

    input_peakdetection = matrix(cbind(x = den$x[in_range], y = den$y[in_range]), ncol = 2)
    colnames(input_peakdetection) = c('x', 'y')

    # Test 5 parametrisations of peakPick neighlim
    pks = Reduce(dplyr::bind_rows,
                 lapply(1:5,
                        function(n) {
                          pk = peakPick::peakpick(mat = input_peakdetection, neighlim = n)
                          input_peakdetection[pk[, 2], , drop = FALSE] %>% as.data.frame()
                        })) %>%
      as_tibble() %>%
      dplyr::arrange(x) %>%
      dplyr::mutate(x = round(x, 2), y = round(y, 2)) %>%
      dplyr::distinct(x, .keep_all = TRUE)

    # library(scorepeak)
    # install.packages('scorepeak')
    #
    #



    hst = hist(mutations$VAF, breaks = seq(0, 1, 0.01), plot = F)$counts
    pks$counts_per_bin = hst[round(pks$x * 100)]

    # Heuristic to remove low-density peaks
    pks = pks %>%
      dplyr::mutate(discarded = y <= max(pks$y) * (1 / 20), from = 'KDE')

    return(list(peaks = pks, density = den))
  }

  # Filter small segments
  candidates = x$snvs$karyotype %>% unique
  candidates = setdiff(candidates,  c("1:1", "1:0", "2:0", "2:1", "2:2", "NA:NA"))
  n_cand = x$n_karyotype[candidates] >= n_min

  analysis = names(n_cand)[n_cand]

  cli::cli_alert_info(
    "Karyotypes {.field {analysis}} with >{.field {n_min}} mutation. Using epsilon = {.field {epsilon}}."
  )

  # Expected peaks
  expected_peaks = lapply(analysis,
                          function(k)
                            expectations_generalised(
                              m = NULL,
                              M = NULL,
                              p = x$purity,
                              karyotype = k
                            )
                          ) %>%
    Reduce(f = bind_rows)

  # Data peaks and densities
  data_fit = lapply(analysis,
                    function(k) {
                      # Get default all data run
                      s_run = aux_peak_calling_general(x$snvs %>% filter(karyotype == k))

                      if (n_bootstrap > 1)
                      {
                        # Bootstrap if required (merge peaks)
                        N = x$snvs %>% filter(karyotype == k) %>% nrow()

                        sb_run = lapply(1:n_bootstrap, function(e) {
                          aux_peak_calling_general(x$snvs %>%
                                                     filter(karyotype == k) %>%
                                                     sample_n(N, replace = TRUE))$peaks
                        }) %>%
                          Reduce(f = bind_rows) %>%
                          distinct(x, .keep_all = TRUE)

                        s_run$peaks = s_run$peaks %>%
                          bind_rows(sb_run) %>%
                          distinct(x, .keep_all = TRUE)
                      }
                      return(s_run)
                    })

  data_densities = lapply(data_fit %>% seq_along,
                          function(f){
                            data.frame(
                                x = data_fit[[f]]$density$x,
                                y = data_fit[[f]]$density$y,
                                karyotype = analysis[f]
                              )
                          }) %>%
    Reduce(f = bind_rows)

  data_peaks = lapply(data_fit %>% seq_along,
                    function(f) {
                      data_fit[[f]]$peaks %>%
                        mutate(karyotype = analysis[f])
                    }) %>%
    Reduce(f = bind_rows)

  # Matching
  for(e in 1:nrow(expected_peaks))
  {
    e_k = expected_peaks$karyotype[e]
    d_k = abs(
      data_peaks %>%
        filter(karyotype == e_k) %>%
        pull(x) -
        expected_peaks$peak[e]
      )

    if(any(d_k < epsilon))
      expected_peaks$matched[e] = TRUE
    else
      expected_peaks$matched[e] = FALSE
  }

  add_counts = function(w) {
    w$n = x$n_karyotype[w$karyotype]
    w %>% as_tibble()
  }

  x$peaks_analysis$general$analysis = analysis
  x$peaks_analysis$general$params = list(n_min = n_min, epsilon = epsilon, n_bootstrap = n_bootstrap)
  x$peaks_analysis$general$expected_peaks = expected_peaks %>% add_counts
  x$peaks_analysis$general$data_peaks = data_peaks %>% add_counts
  x$peaks_analysis$general$data_densities = data_densities %>% add_counts

  x$peaks_analysis$general$summary =
    x$peaks_analysis$general$expected_peaks %>%
      group_by(karyotype, matched, n) %>%
      summarise(hits = n()) %>%
      mutate(matched = ifelse(matched, 'matched', 'mismatched')) %>%
      pivot_wider(names_from = 'matched', values_from = 'hits') %>%
      mutate(matched = ifelse(is.na(matched), 0, matched)) %>%
      mutate(mismatched = ifelse(is.na(mismatched), 0, mismatched)) %>%
      mutate(prop = matched/(matched + mismatched)) %>%
      arrange(prop %>% desc)

  x$peaks_analysis$general$summary %>%
    print()

  return(x)
}

