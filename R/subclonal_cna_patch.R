#' Get clonal peaks
#'
#' @param k simple clonal karyotype
#' @param purity 
#'
#' @return expected peaks
#' @export
#'
#' @examples
get_clonal_peaks <- function(k, purity) {
  multiplicities <- strsplit(k, ":") %>% 
    unlist() %>% 
    as.numeric()
  n_tot <- sum(multiplicities)
  
  if (2 %in% unique(multiplicities)){
    multiplicities <- c(1, 2)}
  
  multiplicities <- multiplicities[multiplicities != 0] 
  
  peaks <- unique(multiplicities) * purity / (n_tot * purity + 2 * (1 - purity))
  return(peaks)
}

#' Get model
#'
#' @param karyotype_combination combination of simple karyotypes of the two subclones
#' @param purity purity to filter result
#' @param ccf 
#'
#' @return dataframe with the possible evolutionary models and corresponding expected peaks
#' @export
#'
#' @examples
get_models = function(karyotype_combination= NULL, purity= NULL, ccf= NULL){
  if (is.null(karyotype_combination)){
    return(CNAqc::model_ccf_purity)
  }else{
    filtered_models = CNAqc::model_ccf_purity[[karyotype_combination]]
    if (!is.null(ccf)){
      lab= paste0(ccf, '-', purity)
      return(filtered_models[[lab]])
    }
  }
}

#' Simple karyotypes 
#'
#' @return a vector of simple karyotypes 
#' @export
#'
#' @examples
simple_karyotypes = function(){
  clonal_karyotypes <- c('1:0', '1:1', '2:0', '2:1', '2:2')
  clonal_karyotypes
}

#' Test whether the expected BAF and DR match the observed ones, i.e. if it is advisable to patch the calls of the segment 
#'
#' @param x a cnaqc object
#' @param seg_id the id of the segment to test
#' @param var expected variance for the BAF distribution
#' @param significance significance level for the Kolmogorov-Smirnov test
#' @param var_dr expected DR rate parameter 
#'
#' @return boolean TRUE or FAlSE, wheter the segment had significantly different BAF of DR distributions
#' @export
#'
#' @examples
test_quality_segment = function(x, seg_id, var = .005, significance=.05, var_dr = 10){
  
  if (length(x$snps %>% filter(segment_id == seg_id) %>% pull(chr))> 5){
    
    cna_call = x$cna %>% filter(segment_id == seg_id)
    n_snps = x$snps %>% filter(segment_id == seg_id) %>% pull(N.BAF) %>% sum()
    # Expected BAF distrbution
    exp_baf <- clonal_expected_baf(paste0(cna_call$Major, ':',cna_call$minor), x$purity)
    if (exp_baf < .01){exp_baf=.01}
    #exp_baf = .01
    #var = .005
    exp_dist_alpha <- ((1 - exp_baf) / var - 1 / exp_baf) * exp_baf ^ 2
    exp_dist_beta <- exp_dist_alpha * (1 / exp_baf - 1)
    e_baf <- rbeta(n_snps, exp_dist_alpha, exp_dist_beta)
    #hist(rbeta(1000, exp_dist_alpha, exp_dist_beta))
    
    # Observed BAF distrubution
    obs_mean_baf = x$snps %>% filter(segment_id == seg_id) %>% pull(BAF) 
    obs_mean_baf = sapply(obs_mean_baf, function(x){
      if(x<.01){
        x=.01
        }
      return(x)
      })
    obs_alpha = ((1 - obs_mean_baf) / var - 1 / obs_mean_baf) * obs_mean_baf ^ 2
    obs_beta = obs_alpha * (1 / obs_mean_baf - 1)
    o_baf = sapply(x$snps %>% filter(segment_id == seg_id) %>% pull(N.BAF), rbeta, shape1=obs_alpha, shape2=obs_beta) %>% unlist()
    
    ks.test(e_baf, o_baf)
    baf_test_p = ks.test(e_baf, o_baf)[["p.value"]]
  
    # Expected DR distribution
    exp_dr <- clonal_expected_dr(paste0(cna_call$Major, ':',cna_call$minor), x$purity)
    exp_dist_shape = exp_dr * var_dr
    e_dr = rgamma(n_snps,exp_dist_shape, var_dr)
    # Observed DR distribution
    obs_mean_dr = x$snps %>% filter(segment_id == seg_id) %>% pull(DR)
    obs_shape = obs_mean_dr * var_dr
    o_dr = sapply(x$snps %>% filter(segment_id == seg_id) %>% pull(N.BAF), rgamma, shape=obs_shape, rate=var_dr) %>% unlist()
    
    ks.test(e_dr, o_dr)
    dr_test_p = ks.test(e_dr, o_dr)[["p.value"]]
  
    if (baf_test_p < significance | dr_test_p < significance){
      TRUE
    }else{
      FALSE
    }
  }else{
    FALSE
  }
  #return(params = list(alpha = alpha, beta = beta))
}

#' Pre-filters the segments with suspect copy number calls based on DR and BAF
#'
#' @param x a cnaqc object
#'
#' @return a vector with the pre-selected segment ids
#' @export
#'
#' @examples
pre_select_segments = function(x){
  segment_ids = x$segment_type %>% filter(segment_type %in% c('simple clonal', 'simple subclonal')) %>% pull(segment_id)
  selected_segs = sapply(segment_ids, test_quality_segment, x=x) 
  selected_segs = selected_segs %>% as.data.frame() %>% mutate(segment_ids=segment_ids)
  colnames(selected_segs) = c('passed', 'segment_ids')
  selected_segs = selected_segs %>% filter(passed == TRUE) %>% pull(segment_ids)
  return(selected_segs)
}

#' Test models for a single segment
#'
#' @param x cnaqc object
#' @param seg_id cnaqc object
#' @param purity 
#' @param top_n how many results to output 
#'
#' @return
#' @export
#'
#' @examples
test_model <- function(x, seg_id, top_n=5, all_solutions=FALSE){
  
  SNV_df = x$mutations %>% filter(segment_id==seg_id)
  SNP_df = x$snps %>% filter(segment_id==seg_id) #%>% filter(!is.na(BAF))
  purity = x$purity
  ploidy = x$ploidy
  
  if (length(SNV_df$chr) ==0 | length(SNP_df$chr) ==0){
    return(NULL)
  }
  
  clonal_results <- clonal_test(SNP_df, SNV_df, purity, ploidy) %>% mutate(segment_id = seg_id) %>% mutate(peak=peaks) %>% select(!peaks)
  sc_results <- sub_clonal_test(SNP_df, SNV_df, purity, ploidy) %>% mutate(segment_id = seg_id)
  
  if (!all_solutions){
      #clonal_results <- clonal_test(SNP_df, SNV_df, purity)
      top_n_clonal = clonal_results %>% mutate(ids = paste0(k1, vaf_ll, baf_ll, dr_ll)) 
      top_n_clonal_ids = top_n_clonal %>% pull(ids) %>% unique()
      top_n_clonal_ids=top_n_clonal_ids[1:top_n]
      top_n_clonal = top_n_clonal %>% filter(ids %in% top_n_clonal_ids) %>% select(!ids)
      clonal_results = top_n_clonal %>% mutate(segment_id = seg_id)
  
      #sc_results <- sub_clonal_test(SNP_df, SNV_df, purity)
      top_n_subclonal = sc_results %>% mutate(ids = paste0(model, ccf_1, vaf_ll)) 
      top_n_subclonal_ids = top_n_subclonal %>% pull(ids) %>% unique() 
      top_n_subclonal_ids =  top_n_subclonal_ids[1:top_n]
      sc_results = top_n_subclonal %>% filter(ids %in% top_n_subclonal_ids) %>% select(!ids) %>% mutate(segment_id = seg_id)
    
    }
  
  seg_cl_baf_ll = clonal_results[1,] %>% pull(baf_ll)
  seg_cl_dr_ll = clonal_results[1,] %>% pull(dr_ll)
  seg_scl_baf_ll = sc_results[1,] %>% pull(baf_ll)
  seg_scl_dr_ll = sc_results[1,] %>% pull(dr_ll)
  
  cl_ll = clonal_results[1,] %>% pull(loglikelihood)
  scl_ll = sc_results[1,] %>% pull(loglikelihood)
  
  if ( ((seg_cl_baf_ll> seg_scl_baf_ll) & (seg_cl_dr_ll> seg_scl_dr_ll)) | (cl_ll > scl_ll) ){
    first_row = clonal_results[1,]
    best = clonal_results %>% filter(k1== first_row$k1, 
                                     k2== first_row$k2, 
                                     vaf_ll== first_row$vaf_ll,
                                     baf_ll==first_row$baf_ll,
                                     dr_ll==first_row$dr_ll, 
                                     loglikelihood==first_row$loglikelihood ) 
  }else{
    first_row = sc_results[1,]
    best = sc_results %>% filter(k1== first_row$k1, 
                                 k2== first_row$k2, 
                                 vaf_ll== first_row$vaf_ll,
                                 baf_ll==first_row$baf_ll,
                                 dr_ll==first_row$dr_ll, 
                                 loglikelihood==first_row$loglikelihood )
  }
  
  #results <- bind_rows(clonal_results, sc_results)
  #results <- results %>% dplyr::arrange(desc(loglikelihood))
  return(list('best'=best, 'clonal'=clonal_results, 'subclonal'=sc_results))
}

#' Clonal test for a segment 
#'
#' @param SNP_df data frame of SNPs sitting on the segment
#' @param SNV_df dataframe of SNVs sitting on the segment
#' @param purity 
#'
#' @return a dataframe of the possible karyotype and corresponding expected peaks given a segment, sorted by their likelihood
#' @export
#'
#' @examples
clonal_test <- function(SNP_df, SNV_df, purity, ploidy=2){
  
  results = lapply(simple_karyotypes(), function(k){
      
    peaks <- data.frame('peaks'=get_clonal_peaks(k, purity))
    # VAF
    v <- VAF_LL(NV= SNV_df$NV, DP= SNV_df$DP, peaks = peaks)
    # BAF
    b = clonal_baf_ll(baf_obs= SNP_df$BAF, n = SNP_df$N.BAF, k1=k, purity=purity)
    # DR
    d <- clonal_dr_ll(SNP_df$DR, n = SNP_df$N.BAF, k, purity, ploidy)
    ll <- v + sum(log(b)) + sum(log(d))
    r <- data.frame(k1 = k, k2 = "", ccf_1 = 1, model = "Clonal", vaf_ll=v, baf_ll=sum(log(b)), dr_ll=sum(log(d)), loglikelihood = ll, peaks=peaks$peaks, model_type= 'Clonal')
    r
  })
  results= Reduce(rbind, results)
  results <- results %>% dplyr::arrange(desc(loglikelihood))
  return(results)
}

#' Subclonal test for a segment 
#'
#' @param SNP_df data frame of SNPs sitting on the segment
#' @param SNV_df dataframe of SNVs sitting on the segment
#' @param purity 
#'
#' @return a dataframe of the possible karyotypes and CCF combinations and corresponding expected peaks given a segment, sorted by their likelihood
#' @export
#'
#' @examples
sub_clonal_test <- function(SNP_df, SNV_df, purity, ploidy=2){
  
  possible_ccf = seq(0.1, 1, by=0.1)
  discrete_ccf = seq(from=0, to= 1, by= .05)
  temp_pur= discrete_ccf[which(abs(discrete_ccf - purity) == min(abs(discrete_ccf - purity)))]
  all_models = names(get_models())
  
  results= lapply(all_models, function(n){
    
    karyotypes <- strsplit(n[1], "-")[[1]]
    k1 <- karyotypes[1]
    k2 <- karyotypes[2]
    
    results_model_n = lapply(possible_ccf, function(ccf){
      df = get_models(karyotype_combination= n, purity= temp_pur, ccf= ccf)
      results_model_nm = data.frame()
      l1 = TRUE
      b1 = TRUE
      for (m in unique(df$model_id)){
        if ('|' %in% strsplit(m, '|')[[1]]){
          if (b1){
            model_type = 'branching 1'
            b1 = FALSE
          }else{
            model_type = 'branching 2'
          }
        }else{
          if (l1){
            model_type = 'linear 1'
            l1 = FALSE
          }else{
            model_type = 'linear 2'
          }
        }
        peak <- df %>% dplyr::filter(model_id == m) %>% dplyr::select(peak) 
        genotypes <- df %>% dplyr::filter(model_id == m) %>% dplyr::select(genotype_1, genotype_2)
        n_peaks <- nrow(peak)
        g1 <- (genotypes$genotype_1 %>% na.omit())[1]
        g2 <- (genotypes$genotype_2 %>% na.omit())[1]
        # VAF likelihood
        v <- VAF_LL(
          NV = SNV_df$NV,
          DP = SNV_df$DP,
          peaks = peak)
        # BAF Likelihood
        b <- BAF_LL(baf_obs = SNP_df$BAF,
                      n = SNP_df$N.BAF,
                      k1 = k1, k2 = k2,
                      purity = purity,
                      ccf = ccf,
                      g1 = g1, 
                      g2 = g2)
        # Depth ratio likelihood
        d <- DR_LL(dpr_obs = SNP_df$DR,
                     n = SNP_df$N.BAF,
                     k1 = k1,
                     k2 = k2,
                     purity = purity, 
                     ccf = ccf,
                   ploidy=ploidy)
        # Overall log-likelihood
        ll <- v + sum(log(b)) + sum(log(d))
        r <- data.frame(k1 = k1, k2 = k2, ccf_1 = ccf, model = m, vaf_ll = v, baf_ll=sum(log(b)), dr_ll=sum(log(d)), loglikelihood = ll, peak=peak, model_type= model_type)
        results_model_nm = rbind(results_model_nm,r)
      }
      results_model_nm
        #results <- bind_rows(results, r)
      })
    results_model_n = Reduce(rbind, results_model_n)
    results_model_n
      
  })
  results= Reduce(rbind, results)
  results <- results %>% dplyr::arrange(desc(loglikelihood))
  return(results)
}


#' Given a vector of segments to test, finds the best copy number solution based on VAF, BAF and DR
#'
#' @param x a cnaqc object
#' @param segments optional, a user-defined vector of segment ids
#' @param top_n optional, alternative to all_solutions, integer declaring how many solutions are to be saved 
#' @param all_solutions optional, alternative to top_n, TRUE/FALSE whether to save all solutions
#' @param preselect optional, whether to pre-select the segment to be tested based on BAF and DR
#'
#' @return the cnaqc object with three new fields : 
#'              - patch_best_solution : the best solutions for each tested segment
#'              - patch_clonal_solution : top n / all clonal solutions for the tested segments
#'              - patch_sunclonal_solution : top n / all subclonal solutions for the tested segments
#' @export
#'
#' @examples
patch = function(x, segments= NULL, top_n=5, all_solutions=FALSE, preselect = FALSE){
  
  if (is.null(segments)){
    segment_ids = x$segment_type %>% filter(segment_type %in% c('simple clonal', 'simple subclonal')) %>% pull(segment_id)
  }else{
    segment_ids = x$segment_type %>% filter(segment_type %in% segments) %>% pull(segment_id)
  }
  
  if (preselect){
    segment_ids = pre_select_segments(x)
  }
  
  cli::cli_h3("Computing solutions for {length(segment_ids)} segments")
  solutions = pbapply::pblapply(segment_ids, test_model, x=x, top_n=top_n, all_solutions=all_solutions)
  solutions <- solutions[!sapply(solutions,is.null)]
  #cli::cli_h3("Segments tesed")
  
  best = lapply(solutions, function(y){
    y$best
  })
  best = Reduce(rbind, best)
  
  x$patch_best_solution = best
  
  clonal_solutions = lapply(solutions, function(y){
    y$clonal
  })
  clonal_solutions = Reduce(rbind, clonal_solutions)
  x$patch_clonal_solutions = clonal_solutions
  
  subclonal_solutions = lapply(solutions, function(y){
    y$subclonal
  })
  subclonal_solutions = Reduce(rbind, subclonal_solutions)
  x$patch_subclonal_solutions = subclonal_solutions
  
  x
  
}

plot_baf_single_segment = function(SNP_df){
  CNAqc:::blank_genome(ref = 'GRCh38', chromosomes = SNP_df$chr[1]) +
    ggplot2::geom_hline(
      yintercept = 0.5,
      color = 'indianred3',
      linetype = 'dashed',
      size = .4
    ) +
    ggplot2::geom_segment(
      data = SNP_df,
      ggplot2::aes(
        x = from,
        xend = to,
        y = BAF,
        yend = BAF
      ),
      size = 1,
      colour = 'black'
    ) +ggplot2::labs(y = 'BAF', x= SNP_df$chr[1])
}
plot_dr_single_segment = function(SNP_df){
  CNAqc:::blank_genome(ref = 'GRCh38', chromosomes = SNP_df$chrs[1]) +
    ggplot2::geom_hline(
      yintercept = 1,
      color = 'indianred3',
      linetype = 'dashed',
      size = .4
    ) +
    ggplot2::geom_segment(
      data = SNP_df %>% filter(DR <= 2.1),
      ggplot2::aes(
        x = from,
        xend = to,
        y = DR,
        yend = DR
      ),
      size = 1,
      colour = 'black'
    ) +
  ggplot2::labs(y = 'Depth ratio', x = SNP_df$chrs[1])
}
plot_vaf_single_segment = function(SNV_df){
  ggplot(SNV_df, aes(x= VAF))+
    geom_histogram(bins=100, fill= 'forestgreen', alpha= .5) +
    my_ggplot_theme() +
    labs(x='VAF',y='')
}
plot_patch_best_solution = function(x, seg_id){
  snvs_seg = x$mutations %>% filter(segment_id == seg_id)
  snps_seg = x$snps %>% filter(segment_id == seg_id)
  cna_seg = x$cna %>% filter(segment_id == seg_id)
  
  best_solution = x$patch_best_solution %>% filter(segment_id == seg_id)
  
  e_baf = clonal_expected_baf(paste0(cna_seg$Major, ':', cna_seg$minor), x$purity)
  e_dr = clonal_expected_dr(paste0(cna_seg$Major, ':', cna_seg$minor), x$purity, x$ploidy)
  e_peaks = get_clonal_peaks(paste0(cna_seg$Major, ':', cna_seg$minor), x$purity)
  
  if (x$patch_best_solution$model[1] == 'Clonal'){
    
    proposed_baf = clonal_expected_baf(best_solution$k1, x$purity)
    proposed_dr = clonal_expected_dr(best_solution$k1, x$purity, x$ploidy)
    proposed_peaks = get_clonal_peaks(best_solution$k1, x$purity)
    
  }else{
    k1 = best_solution$k1 %>% unique()
    k2 = best_solution$k2 %>% unique()
    ccf = best_solution$ccf_1 %>% unique()
    g1 = strsplit(best_solution$model[1], '->')[[1]][1]
    g2 = strsplit(best_solution$model[1], '->')[[1]][2]
    
    proposed_baf = expected_baf(k1, k2, x$purity, ccf, g1, g2)
    proposed_dr = expected_dr(k1, k2, x$purity, ccf, x$ploidy)
    proposed_peaks = best_solution$peak
  }
  
  
  old_baf_plot = plot_baf_single_segment(snps_seg) + geom_hline(yintercept = e_baf, color= 'forestgreen')
  new_baf_plot = plot_baf_single_segment(snps_seg) + geom_hline(yintercept = proposed_baf, color= 'forestgreen')
  
  old_dr_plot = plot_dr_single_segment(snps_seg) + geom_hline(yintercept = e_dr, color= 'forestgreen')
  new_dr_plot = plot_dr_single_segment(snps_seg) + geom_hline(yintercept = proposed_dr, color= 'forestgreen')
  
  old_vaf = plot_vaf_single_segment(snvs_seg) + geom_vline(xintercept = e_peaks, color = 'forestgreen', linetype= 'dashed')
  new_vaf = plot_vaf_single_segment(snvs_seg) + geom_vline(xintercept = proposed_peaks, color = 'forestgreen', linetype= 'dashed')
  
  st= 'AAAACC
       BBBBCC'
  old_solution_name = paste0('Original solution: clonal segment of karyotype ', paste0(cna_seg$Major, ':', cna_seg$minor))
  old_plot = patchwork::wrap_plots(old_baf_plot, old_dr_plot, old_vaf, design=st) +  patchwork::plot_annotation(title = old_solution_name) & 
    theme(text = element_text(size = 10))
    #ggplot2::ggtitle(old_solution_name)
  which_solution = 'subclonal'
  if (x$patch_best_solution$model[1] == 'Clonal'){
    which_solution = 'clonal'
  }
  new_solution_name = paste0('Proposed solution: ', which_solution, ' segment with karyotype ', best_solution$k1 %>% unique() , ' ', best_solution$k2 %>% unique(), ', CCF ', ccf = best_solution$ccf_1 %>% unique())          
  new_plot = patchwork::wrap_plots(new_baf_plot, new_dr_plot, new_vaf, design=st) +  patchwork::plot_annotation(title = new_solution_name) &  theme(text = element_text(size = 10))
  
  final_plot = ggarrange(plotlist = list(old_plot, new_plot), ncol= 1)
}

plot_patch_all_solutions_subclonal = function(x, seg_id){
  snvs_seg = x$mutations %>% filter(segment_id == seg_id)
  snps_seg = x$snps %>% filter(segment_id == seg_id)
  cna_seg = x$cna %>% filter(segment_id == seg_id)
  
  subclonal_solutions = x$patch_subclonal_solutions %>% filter(segment_id == seg_id)
  
  subclonal_solutions_plot = subclonal_solutions %>% 
    ggplot(aes(x = ccf_1, y = loglikelihood, color= model_type)) +
    geom_line() +
    scale_color_manual(values = c('branching 1' = '#D741A7', 'branching 2' = '#3A1772', 'linear 1'= '#F2CD5D', 'linear 2'= '#DEA54B'))+
    ggh4x::facet_grid2(k1~k2, scales = 'free') + my_ggplot_theme() + labs(x= 'CCF') +
    guides(color = guide_legend(title = "Model type"))
  
  library(grid)
  library(gtable)
  grob <- ggplot2::ggplotGrob(subclonal_solutions_plot)
  idx <- which(grob$layout$name %in% c("panel-1-2", "panel-1-3", "panel-1-4", "panel-1-5",
                                       "panel-2-3", "panel-2-4", "panel-2-5",
                                       "panel-3-4", "panel-3-5", 
                                       "panel-4-5"))
  for (i in idx) grob$grobs[[i]] <- nullGrob()
  
  # Move x axes up
  idx <- which(grob$layout$name %in% c("axis-b-1-5"))
  grob$layout[idx, c("t", "b")] <- grob$layout[idx, c("t", "b")] - c(14, 2)
  
  idx <- which(grob$layout$name %in% c("axis-b-2-5"))
  grob$layout[idx, c("t", "b")] <- grob$layout[idx, c("t", "b")] - c(12, 2)
  
  idx <- which(grob$layout$name %in% c("axis-b-3-5"))
  grob$layout[idx, c("t", "b")] <- grob$layout[idx, c("t", "b")] - c(8, 2)
  
  idx <- which(grob$layout$name %in% c("axis-b-4-5"))
  grob$layout[idx, c("t", "b")] <- grob$layout[idx, c("t", "b")] - c(4, 2)
  
  # Move y axes left
  idx <- which(grob$layout$name %in% c("axis-l-2-1"))
  grob$layout[idx, c("l", "r")] <- grob$layout[idx, c("l", "r")] + c(2, 4)
  
  idx <- which(grob$layout$name %in% c("axis-l-3-1"))
  grob$layout[idx, c("l", "r")] <- grob$layout[idx, c("l", "r")] + c(2, 8)
  
  idx <- which(grob$layout$name %in% c("axis-l-4-1"))
  grob$layout[idx, c("l", "r")] <- grob$layout[idx, c("l", "r")] + c(2, 12)
  
  idx <- which(grob$layout$name %in% c("axis-l-5-1"))
  grob$layout[idx, c("l", "r")] <- grob$layout[idx, c("l", "r")] + c(2, 16)
  
  #grid.newpage()
  #grid.draw(grob)
  ggpubr::as_ggplot(grob) + ggtitle('Subclonal solutions')
}



patch_plot = function(x,seg_id){
  ggarrange(plotlist= list(plot_patch_best_solution(x, seg_id), plot_patch_all_solutions_subclonal(x,seg_id), ncol=2))
}




