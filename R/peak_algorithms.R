###### ###### ############ ###### ############ ###### ############ ###### ######
# SIMPLE KARYOTYPES
analyze_peaks_common = function(x,
                                karyotypes = c('1:0', '1:1', '2:0', '2:1', '2:2'),
                                min_karyotype_size = 0,
                                min_absolute_karyotype_mutations = 100,
                                p_binsize_peaks = 0.005,
                                purity_error = 0.05,
                                VAF_tolerance = 0.015,
                                n_bootstrap = 1,
                                kernel_adjust = 1,
                                matching_strategy = "closest",
                                KDE = TRUE, 
                                min_VAF = 0)
{
  # Karyotypes of interest, and filter for karyotype size
  analysis_type = x$n_karyotype[names(x$n_karyotype) %in%  karyotypes] %>% names()
  analysis_numb = x$n_karyotype[x$n_karyotype >= min_absolute_karyotype_mutations] %>% names()
  analysis_prop = x$n_karyotype[x$n_karyotype/sum(x$n_karyotype) >= min_karyotype_size] %>% names()
  
  analysis = intersect(analysis_type, analysis_numb) %>% intersect(analysis_prop)
  
  if (length(analysis) == 0) {
    cli::cli_alert_warning("No karyotypes satisfy input data filters.")
    return(x)
  }
  else
  {
    n = x$n_karyotype[names(x$n_karyotype) %in% analysis] %>% sum()
    
    cli::cli_alert_info(
      paste0(
        "Analysing {.field {n}} mutations mapping to karyotype(s) {.field {analysis}}."
      )
    )
  }
  
  # Expected peaks
  expected_peaks = lapply(analysis,
                          function(k)
                          {
                            AB = as.numeric(strsplit(k, ':')[[1]])
                            CNAqc:::expected_vaf_peak(AB[1], AB[2], x$purity)
                          }) %>%
    Reduce(f = bind_rows) %>%
    left_join(
      delta_vaf_karyo(epsilon_error = purity_error,
                      purity = x$purity) %>%
        filter(karyotype %in% analysis),
      by = c('karyotype', 'mutation_multiplicity')
    )
  
  # Run peak detection
  data_fits = x$mutations %>%
    dplyr::filter(VAF > min_VAF) %>% 
    dplyr::filter(karyotype %in% analysis) %>%
    dplyr::group_split(karyotype) %>%
    lapply(
      FUN = function(w) {
        cli::cli_alert_info("Mixed type peak detection for karyotype {.field {w$karyotype[1]}} ({.field {x$n_karyotype[w$karyotype[1]]}} mutations)")
        
        w %>% combined_peak_detector(kernel_adjust = kernel_adjust,
                                     n_bootstrap = n_bootstrap)
      }
    )
  
  names(data_fits) = x$mutations %>%
    dplyr::filter(VAF > min_VAF) %>% 
    dplyr::filter(karyotype %in% analysis) %>%
    dplyr::group_split(karyotype) %>%
    sapply(function(e) e$karyotype[1])
  
  # Weights by karyotype
  analysis_weight = x$n_karyotype[names(x$n_karyotype) %in% analysis]
  analysis_weight = analysis_weight/sum(analysis_weight)
  
  # Match peaks
  control = function(k)
  {
    peaks = data_fits[[k]]
    expectation = expected_peaks %>% filter(karyotype == k)
    
    # =-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Compute matching ~ get any possible match given the expectation
    #
    # 1) sort the expected peaks and the data peaks
    # 2) combine them and subtract
    # =-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    get_match = function(p, peaks)
    {
      not_discarded = peaks %>% filter(!discarded)
      
      distances = abs(not_discarded$x - expectation$peak[p])
      id_match = which.min(distances)
      
      if (length(id_match) == 0)
        return(NA)
      else
        return(not_discarded[id_match,])
    }
    
    matched_peaks = lapply(seq_along(expectation$peak), get_match, peaks = peaks$peaks)
    
    # Matching table
    matching = expectation %>%
      dplyr::bind_cols(Reduce(dplyr::bind_rows, matched_peaks))
    
    # Distance in VAF space, converted to purity space, and matched with bands
    matching = matching %>%
      rowwise() %>%
      mutate(
        offset_VAF = peak - x,
        # VAF space
        offset = compute_delta_purity(
          vaf = x,
          delta_vaf = offset_VAF,
          ploidy = strsplit(k, ':') %>% unlist %>% as.numeric() %>% sum(),
          multiplicity = mutation_multiplicity
        ),
        weight = analysis_weight[k],
        epsilon = purity_error,
        VAF_tolerance = VAF_tolerance,
        matched = overlap_bands(
          peak = x,
          tolerance = VAF_tolerance,
          left_extremum = peak - delta_vaf,
          right_extremum = peak + delta_vaf
        )
      )
  }
  
  qc_table = lapply(analysis, control) %>% Reduce(f = bind_rows)
  
  # =-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # Compute a linear combination for a GOF of the sample and each karyotype
  # =-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  overall_score = qc_table$weight %*% qc_table$offset %>% as.numeric()
  karyotype_score = qc_table %>% group_by(karyotype) %>% summarise(score = sum(weight * offset))
  
  # QC per karyotype
  qc_per_karyotype = qc_table %>%
    dplyr::group_by(karyotype) %>%
    dplyr::arrange(desc(counts_per_bin)) %>%
    dplyr::filter(row_number() == 1) %>%
    dplyr::select(karyotype, matched) %>%
    dplyr::mutate(QC = ifelse(matched, "PASS", 'FAIL')) %>%
    dplyr::select(-matched)
  
  qc_table = qc_table %>% dplyr::left_join(qc_per_karyotype, by = 'karyotype')
  
  print(qc_table)
  
  # Results assembly
  QC = qc_table %>%
    dplyr::group_by(QC) %>%
    dplyr::summarise(prop = sum(weight)) %>%
    dplyr::arrange(desc(prop)) %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::pull(QC)
  
  if (QC == "FAIL")
    cli::cli_alert_danger(
      "Peak detection {red('FAIL')} with {.value {red(paste0('r = ', overall_score))}} - maximum purity error \u03B5 = {.field {purity_error}}."
    )
  
  if (QC == "PASS")
    cli::cli_alert_success(
      "Peak detection {green('PASS')} with {.value {green(paste0('r = ', overall_score))}} - maximum purity error \u03B5 = {.field {purity_error}}."
    )
  
  fits = names(data_fits) %>%
    lapply(function(k){
      list(
        xy_peaks = data_fits[[k]]$peaks,
        density = data_fits[[k]]$density,
        matching = qc_table %>% filter(karyotype == k)
      )
    })
  names(fits) = names(data_fits)
  
  x$peaks_analysis = list(
    score = overall_score,
    fits = fits,
    matches = qc_table,
    matching_strategy = matching_strategy,
    purity_error = purity_error,
    # plots = plots,
    min_karyotype_size = min_karyotype_size,
    p_binsize_peaks = p_binsize_peaks,
    # matching_epsilon = matching_epsilon,
    QC = QC,
    KDE = NA
  )
  
  return(x)
}
###### ###### ############ ###### ############ ###### ############ ###### ######

###### ###### ############ ###### ############ ###### ############ ###### ######
# RARE COMPLEX KARYOTYPES
analyze_peaks_general = function(x,
                                 n_min = 50,
                                 epsilon = 0.03,
                                 kernel_adjust = 1,
                                 n_bootstrap = 5, 
                                 min_VAF = 0)
{
  # Filter small segments
  candidates = x$mutations%>% dplyr::filter(VAF > min_VAF) %>% dplyr::pull(karyotype) %>% unique
  candidates = setdiff(candidates,  c("1:1", "1:0", "2:0", "2:1", "2:2", "NA:NA"))
  n_cand = x$n_karyotype[candidates] >= n_min
  
  analysis = names(n_cand)[n_cand]
  
  # Expected peaks
  expected_peaks = lapply(analysis,
                          function(k)
                            expectations_generalised(
                              m = NULL,
                              M = NULL,
                              p = x$purity,
                              karyotype = k
                            )) %>%
    Reduce(f = bind_rows)
  
  # Data peaks and densities
  data_fit = x$mutations %>%
    dplyr::filter(VAF > min_VAF) %>% 
    filter(karyotype %in% analysis) %>%
    group_split(karyotype) %>%
    lapply(simple_peak_detector,
           kernel_adjust = kernel_adjust,
           n_bootstrap = n_bootstrap)
  
  names(data_fit) = x$mutations %>%
    dplyr::filter(VAF > min_VAF) %>% 
    filter(karyotype %in% analysis) %>%
    group_split(karyotype) %>%
    sapply(function(e) e$karyotype[1])
  
  
  data_densities = lapply(data_fit %>% names,
                          function(f) {
                            data.frame(
                              x = data_fit[[f]]$density$x,
                              y = data_fit[[f]]$density$y,
                              karyotype = f
                            )
                          }) %>%
    Reduce(f = bind_rows)
  
  data_peaks = lapply(data_fit %>% names,
                      function(f) {
                        data_fit[[f]]$peaks %>%
                          mutate(karyotype = f)
                      }) %>%
    Reduce(f = bind_rows)
  
  # Matching
  for (e in 1:nrow(expected_peaks))
  {
    e_k = expected_peaks$karyotype[e]
    d_k = abs(data_peaks %>%
                filter(karyotype == e_k) %>%
                pull(x) -
                expected_peaks$peak[e])
    
    if (any(d_k < epsilon))
      expected_peaks$matched[e] = TRUE
    else
      expected_peaks$matched[e] = FALSE
  }
  
  add_counts = function(w) {
    w$n = x$n_karyotype[w$karyotype]
    w %>% as_tibble()
  }
  
  # Results
  x$peaks_analysis$general$analysis = analysis
  x$peaks_analysis$general$params = list(n_min = n_min,
                                         epsilon = epsilon,
                                         n_bootstrap = n_bootstrap)
  x$peaks_analysis$general$expected_peaks = expected_peaks %>% add_counts
  x$peaks_analysis$general$data_peaks = data_peaks %>% add_counts
  x$peaks_analysis$general$data_densities = data_densities %>% add_counts
  
  # Summary table
  stable =
    x$peaks_analysis$general$expected_peaks %>%
    group_by(karyotype, matched, n) %>%
    summarise(hits = n()) %>%
    mutate(matched = ifelse(matched, 'matched', 'mismatched')) %>%
    pivot_wider(names_from = 'matched', values_from = 'hits')
  
  if ("matched" %in% (stable %>% colnames))
    stable = stable %>%
    mutate(matched = ifelse(is.na(matched), 0, matched))
  else
    stable$matched = 0
  
  if ("mismatched" %in% (stable %>% colnames))
    stable = stable %>%
    mutate(mismatched = ifelse(is.na(mismatched), 0, mismatched))
  else
    stable$mismatched = 0
  
  stable = stable %>%
    mutate(prop = matched / (matched + mismatched)) %>%
    arrange(dplyr::desc(prop))
  
  x$peaks_analysis$general$summary = stable
  
  return(x)
}
###### ###### ############ ###### ############ ###### ############ ###### ######


###### ###### ############ ###### ############ ###### ############ ###### ######
# SUBCLONAL PEAKS
analyze_peaks_subclonal = function(x,
                                   epsilon = 0.025,
                                   n_min = 100,
                                   n_bootstrap = 5,
                                   kernel_adjust = 1,
                                   starting_state = '1:1',
                                   cluster_subclonal_CCF = FALSE, 
                                   min_VAF = 0)
{
  if (is.null(x$cna_subclonal) | nrow(x$cna_subclonal) == 0)
    return(x)
  
  # What we can inspect
  subclonal_calls = x$cna_subclonal %>% filter(n > n_min)
  
  if (nrow(subclonal_calls) == 0)
    return(x)
  
  # if we cluster CCF in subclones than we don't run sample by sample
  if(cluster_subclonal_CCF & nrow(subclonal_calls) > 1) {
    
    # cluster mutations
    clusts <- Mclust(subclonal_calls$CCF,1:20, modelNames = "E")
    
    cli::cli_alert(
      "Clustered CCFs into {.field {clusts$G}} subclones with BIC={.field {round(clusts$bic,2)}}"
    )
    
    subclonal_calls$mclusters <- clusts$classification
    # aggregate mutations with the same clusters
    subclonal_calls_old <- subclonal_calls
    
    subclonal_calls <- subclonal_calls %>% group_by(mclusters, karyotype,karyotype_2) %>% 
      summarize(CCF = mean(CCF),
                n = sum(n), mutations = list(do.call(rbind, mutations))) 
    
    # Add segment id to mutations data
    subclonal_mutations = lapply(1:nrow(subclonal_calls),
                                 function(i)
                                 {
                                   CCF_1 = subclonal_calls$CCF[i] %>% round(2)
                                   karyotype_1 = subclonal_calls$karyotype[i]
                                   karyotype_2 = subclonal_calls$karyotype_2[i]
                                   mclusters = subclonal_calls$mclusters[i]
                                   n = subclonal_calls$n[i]
                                   
                                   subclonal_calls$mutations[[i]] %>% 
                                   dplyr::filter(VAF > min_VAF) %>% 
                                   mutate(
                                     segment_id = paste0(
                                       'Subclone ', mclusters,
                                       " ( ",karyotype_1, " | ", karyotype_2,
                                       " ) \n n = ",
                                       n,
                                       ' \n',
                                       karyotype_1,
                                       ' ',
                                       CCF_1,
                                       ' ',
                                       karyotype_2,
                                       ' ',
                                       1 - CCF_1)
                                   )
                                 }) %>%
      Reduce(f = bind_rows)
    
    # Our QC atom is the segment 
  } else {
    
    # Add segment id to mutations data
    subclonal_mutations = lapply(1:nrow(subclonal_calls),
                                 function(i)
                                 {
                                   CCF_1 = subclonal_calls$CCF[i] %>% round(2)
                                   karyotype_1 = subclonal_calls$karyotype[i]
                                   karyotype_2 = subclonal_calls$karyotype_2[i]
                                   
                                   L = round((subclonal_calls$to[i] - subclonal_calls$from[i]) /
                                               1e6, 1)
                                   
                                   subclonal_calls$mutations[[i]] %>%  
                                    dplyr::filter(VAF > min_VAF) %>% 
                                    mutate(
                                     segment_id = paste0(
                                       subclonal_calls$chr[i],
                                       ' ',
                                       subclonal_calls$from[i],
                                       ' ',
                                       "\n(",
                                       L,
                                       "Mb, n = ",
                                       subclonal_calls$n[i],
                                       ')\n',
                                       karyotype_1,
                                       ' ',
                                       CCF_1,
                                       ' ',
                                       karyotype_2,
                                       ' ',
                                       1 - CCF_1
                                     )
                                   )
                                 }) %>%
      Reduce(f = bind_rows)
    
    
  }
  
  
  
  all_segments = subclonal_mutations$segment_id %>% unique()
  
  cli::cli_alert(
    "Computing evolution models for subclonal CNAs - starting from {.field {starting_state}}"
  )
  
  expected_peaks = easypar::run(
    FUN = function(i) {
      expectations_subclonal(
        starting = starting_state,
        CCF_1 = subclonal_calls$CCF[i],
        karyotype_1 = subclonal_calls$karyotype[i],
        karyotype_2 = subclonal_calls$karyotype_2[i],
        purity = x$purity
      ) %>%
        mutate(segment_id = (subclonal_mutations$segment_id %>% unique)[i])
      
    },
    PARAMS = lapply(1:nrow(subclonal_calls), list),
    parallel = FALSE
  )
  
  expected_peaks = Reduce(bind_rows, expected_peaks)
  
  all(expected_peaks$segment_id %>% unique() %in% all_segments)
  
  # Data peaks and densities
  data_fit = easypar::run(
    FUN = function(i) {
      subclonal_calls$mutations[[i]] %>%
        simple_peak_detector(kernel_adjust = kernel_adjust, n_bootstrap = n_bootstrap)
    },
    PARAMS = lapply(1:nrow(subclonal_calls), list),
    parallel = FALSE,
    filter_errors = FALSE
  )
  
  # Densities
  data_densities = lapply(data_fit %>% seq_along,
                          function(f) {
                            data.frame(
                              x = data_fit[[f]]$density$x,
                              y = data_fit[[f]]$density$y,
                              segment_id = (subclonal_mutations$segment_id %>% unique)[f]
                            )
                          }) %>%
    Reduce(f = bind_rows)
  
  # Peaks
  data_peaks = lapply(data_fit %>% seq_along,
                      function(f) {
                        data_fit[[f]]$peaks %>%
                          mutate(segment_id = (subclonal_mutations$segment_id %>% unique)[f])
                      }) %>%
    Reduce(f = bind_rows)
  
  # Matching
  expected_peaks$matched = sapply(1:nrow(expected_peaks), function(i)
  {
    t_i = expected_peaks$peak[i]
    d_i = data_peaks %>%
      filter(segment_id == expected_peaks$segment_id[i])
    
    return(any(abs(d_i$x - t_i) <= epsilon))
  })
  
  # Decide the mode of evolution if possible
  decision_table = NULL
  
  for (s in expected_peaks$segment_id %>% unique) {
    rankings = expected_peaks %>%
      filter(segment_id == s) %>%
      group_by(model_id) %>%
      summarise(prop = sum(matched == TRUE) / n()) %>%
      arrange(desc(prop))
    
    best_choice = rankings %>% filter(prop == rankings$prop[1]) %>% pull(model_id)
    
    rankings = rankings %>%
      filter(model_id %in% best_choice) %>%
      mutate(segment_id = s,
             model = ifelse(grepl("\\|", model_id), "branching", "linear")) %>%
      select(segment_id, model_id, model, prop)
    
    decision_table = bind_rows(decision_table, rankings)
  }
  
  if(cluster_subclonal_CCF & nrow(subclonal_calls) > 1){
    splitter = function(x)
    {
      x$size = strsplit(x$segment_id, split = '\\n') %>%
        sapply(function(x)
          x[[2]])
      x$clones = strsplit(x$segment_id, split = '\\n') %>%
        sapply(function(x)
          x[[3]])
      x$segment_id = strsplit(x$segment_id, split = '\\n') %>%
        sapply(function(x)
          x[[1]]) 
      
      x
    }
  } else {
    splitter = function(x)
    {
      x$size = strsplit(x$segment_id, split = '\\n') %>%
        sapply(function(x)
          x[[2]])
      x$clones = strsplit(x$segment_id, split = '\\n') %>%
        sapply(function(x)
          x[[3]])
      x$segment_id = strsplit(x$segment_id, split = '\\n') %>%
        sapply(function(x)
          x[[1]]) %>%
        strsplit(split = ' ') %>%
        sapply(function(x) {
          paste(x[1], x[2], sep = '@')
        })
      
      x
    }
  }
  
  
  # Results
  x$peaks_analysis$subclonal$params = list(n_min = n_min,
                                           epsilon = epsilon,
                                           n_bootstrap = n_bootstrap,
                                           clustered = cluster_subclonal_CCF)
  x$peaks_analysis$subclonal$expected_peaks = expected_peaks %>% splitter
  x$peaks_analysis$subclonal$data_peaks = data_peaks %>% splitter
  x$peaks_analysis$subclonal$data_densities = data_densities %>% splitter
  x$peaks_analysis$subclonal$mutations = subclonal_mutations %>% splitter
  x$peaks_analysis$subclonal$summary = decision_table %>% splitter
  
  if(cluster_subclonal_CCF & nrow(subclonal_calls) > 1) 
    x$peaks_analysis$subclonal$mclust <- clusts
  
  return(x)
}
###### ###### ############ ###### ############ ###### ############ ###### ######


###### ###### ############ ###### ############ ###### ############ ###### ######
# Peak matching functions - density smoothing
smooth_data = function(mutations, kernel_adjust)
{
  # Smoothed Gaussian kernel for VAF
  y = mutations %>% dplyr::pull(VAF)
  
  density(y,
          kernel = 'gaussian',
          adjust = kernel_adjust,
          na.rm = T)
}
###### ###### ############ ###### ############ ###### ############ ###### ######

# Peak matching functions - phase against density smoothing
phase_to_density = function(peaks, density)
{
  # Adjust peaks height based on KDE
  target_density = density$x
  
  peaks %>%
    dplyr::rowwise() %>%
    dplyr::mutate(which_x = which.min(abs(target_density - x)),
                  y = density$y[which_x]) %>%
    dplyr::select(-which_x)
}
###### ###### ############ ###### ############ ###### ############ ###### ######

# Peak matching functions - match by KDE
simple_peak_detector = function(mutations, kernel_adjust, n_bootstrap)
{
  # Single run
  single_run = function(mutations, ...)
  {
    xy_peaks = den = NULL
    
    # Smoothed Gaussian kernel for VAF
    den = mutations %>% smooth_data(...)
    
    # in_range = den$x >= min(y, na.rm = T) & den$x <= max(y, na.rm = T)
    in_range = TRUE
    
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
      dplyr::distinct(x, .keep_all = TRUE) %>%
      dplyr::mutate(x = case_when(x > 1 & x < 1.01 ~ 1,
                                  x < 0 & x > -0.01 ~ 0,
                                  TRUE ~ x),) %>%
      dplyr::filter(x <= 1, x >= 0)
    
    hst = hist(mutations$VAF,
               breaks = seq(0, 1, 0.01),
               plot = F)$counts
    
    indexes = round(pks$x * 100)
    if (indexes[1]==0){indexes[1]=1}
    pks$counts_per_bin = hst[indexes]
    
    # pks$counts_per_bin = hst[round(pks$x * 100)]
    
    # Heuristic to remove low-density peaks
    pks = pks %>%
      dplyr::mutate(discarded = y <= max(pks$y) * (1 / 20),
                    from = 'KDE')
    
    return(list(peaks = pks, density = den))
    
  }
  
  # Get default: all data run
  s_run = mutations %>% single_run(kernel_adjust)
  
  # Bootstrap and merge peaks
  if (n_bootstrap > 1)
  {
    sb_run = lapply(1:n_bootstrap, function(e) {
      single_run(mutations %>% sample_n(mutations %>% nrow(), replace = TRUE),
                 kernel_adjust = kernel_adjust)$peaks
    }) %>%
      Reduce(f = bind_rows) %>%
      distinct(x, .keep_all = TRUE) %>%
      phase_to_density(density = s_run$density)
    
    s_run$peaks = s_run$peaks %>%
      bind_rows(sb_run) %>%
      distinct(x, .keep_all = TRUE)
  }
  
  return(s_run)
}
###### ###### ############ ###### ############ ###### ############ ###### ######

# Peak matching functions - match by mixture model
mixture_peak_detector = function(mutations, kernel_adjust, n_bootstrap)
{
  if (mutations$VAF %>% unique() %>% length() == 1)
    return(NULL)
  
  single_run = function(mutations, ...)
  {
    # BMix clustering
    bm = BMix::bmixfit(
      data.frame(successes = mutations$NV, trials = mutations$DP),
      K.BetaBinomials = 0,
      K.Binomials = 1:4,
      silent = TRUE
    )
    
    # Smoothed Gaussian kernel for VAF
    den = mutations %>% smooth_data(...)
    hst = hist(mutations$VAF,
               breaks = seq(0, 1, 0.01),
               plot = F)$counts
    
    llxy = NULL
    for (b in names(bm$B.params))
    {
      w_den = which.min(abs(den$x - bm$B.params[b]))
      
      tnw = tibble(x = den$x[w_den],
                   y = den$y[w_den],
                   # counts_per_bin = bm$pi[b] * (mutations %>% nrow), # Wrong
                   discarded = FALSE)
      
      # Counts are counted the same way regardless it is a BMix fit or not.
      tnw$counts_per_bin = hst[round(tnw$x * 100)]
      
      llxy = llxy %>%
        bind_rows(tnw)
    }
    
    llxy %>% mutate(from = 'BMix') %>% return()
  }
  
  single_run_kde = function(mutations, ...)
  {
    xy_peaks = den = NULL
    
    # Smoothed Gaussian kernel for VAF
    den = mutations %>% smooth_data(...)
    
    # in_range = den$x >= min(y, na.rm = T) & den$x <= max(y, na.rm = T)
    in_range = TRUE
    
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
      dplyr::distinct(x, .keep_all = TRUE) %>%
      dplyr::mutate(x = case_when(x > 1 & x < 1.01 ~ 1,
                                  x < 0 & x > -0.01 ~ 0,
                                  TRUE ~ x),) %>%
      dplyr::filter(x <= 1, x >= 0)
    
    hst = hist(mutations$VAF,
               breaks = seq(0, 1, 0.01),
               plot = F)$counts
    pks$counts_per_bin = hst[round(pks$x * 100)]
    
    # Heuristic to remove low-density peaks
    pks = pks %>%
      dplyr::mutate(discarded = y <= max(pks$y) * (1 / 20),
                    from = 'KDE')
    
    return(list(peaks = pks, density = den))
    
  }
  
  # Get default: all data run
  s_run = mutations %>% single_run(kernel_adjust)
  
  # Bootstrap and merge peaks
  if (n_bootstrap > 1)
  {
    sb_run = lapply(1:n_bootstrap, function(e) {
      single_run(mutations %>% sample_n(mutations %>% nrow(), replace = TRUE),
                 kernel_adjust = kernel_adjust)
    }) %>%
      Reduce(f = bind_rows) %>%
      distinct(x, .keep_all = TRUE) 
    
    s_kde = single_run_kde(mutations %>% sample_n(mutations %>% nrow(), replace = TRUE),
                           kernel_adjust = kernel_adjust)
    s_run = sb_run %>%
      phase_to_density(density = s_kde$density)
    
    s_run = s_run %>%
      #bind_rows(sb_run) %>%
      distinct(x, .keep_all = TRUE)
  }
  
  return(s_run)
}
###### ###### ############ ###### ############ ###### ############ ###### ######

# Peak matching functions - match by KDE + mixture model
combined_peak_detector = function(mutations, kernel_adjust, n_bootstrap)
{
  kde_peaks = mutations %>%
    simple_peak_detector(kernel_adjust = kernel_adjust, n_bootstrap = n_bootstrap)
  
  mixture_peaks = mutations %>%
    mixture_peak_detector(kernel_adjust = kernel_adjust, n_bootstrap = n_bootstrap)
  
  return(
    list(
      peaks = kde_peaks$peaks %>% bind_rows(mixture_peaks),
      density = kde_peaks$density
    )
  )
}
###### ###### ############ ###### ############ ###### ############ ###### ######