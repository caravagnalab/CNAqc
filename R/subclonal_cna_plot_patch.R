#' Plot the BAF of bins sitting on a selected segment
#'
#' @param SNP_df dataframe with the binned SNPs sitting on the segment 
#'
#' @return 
#' @export
#'
#' @examples
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

#' Plot the DR of bins sitting on a selected segment
#'
#' @param SNP_df dataframe with the binned SNPs sitting on the segment 
#'
#' @return
#' @export
#'
#' @examples
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


#' Plot the VAF of SNVs sitting on a single segment
#'
#' @param SNV_df dataframe of SNVs sitting on the segment 
#'
#' @return
#' @export
#'
#' @examples
plot_vaf_single_segment = function(SNV_df){
  ggplot(SNV_df, aes(x= VAF))+
    geom_histogram(bins=100, fill= 'forestgreen', alpha= .5) +
    my_ggplot_theme() +
    labs(x='VAF',y='')
}

#' Plot comparing the original solution with the best one found by function patch (BAF, DR, VAF)
#'
#' @param x cnaqc object
#' @param seg_id id of the segment of interest
#'
#' @return
#' @export
#'
#' @examples
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
  
  st= 'AAACC
       BBBCC'
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

#' Plot summarising all subclonal solutions for a single segment found by patch in terms of loglikelihood and CCF, dividing them by karyotype combination and model type
#'
#' @param x a cnaqc object
#' @param seg_id the segment of interest
#'
#' @return
#' @export
#'
#' @examples
plot_patch_all_solutions_subclonal = function(x, seg_id){
  snvs_seg = x$mutations %>% filter(segment_id == seg_id)
  snps_seg = x$snps %>% filter(segment_id == seg_id)
  cna_seg = x$cna %>% filter(segment_id == seg_id)
  
  subclonal_solutions = x$patch_subclonal_solutions %>% filter(segment_id == seg_id)
  
  subclonal_solutions_plot = subclonal_solutions %>% 
    ggplot(aes(x = ccf_1, y = loglikelihood, color= model_type)) +
    geom_point(size=1) +
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
  ggpubr::as_ggplot(grob) + ggplot2::ggtitle('Subclonal solutions')
}

#' Plot summarising the clonal solutions for a single segment found by patch 
#'
#' @param x a cnaqc object
#' @param seg_id 
#'
#' @return
#' @export
#'
#' @examples
plot_patch_all_solutions_clonal = function(x, seg_id){
  snvs_seg = x$mutations %>% filter(segment_id == seg_id)
  snps_seg = x$snps %>% filter(segment_id == seg_id)
  cna_seg = x$cna %>% filter(segment_id == seg_id)
  
  clonal_solutions = x$patch_clonal_solutions %>% filter(segment_id == seg_id)
  
  clonal_solutions_plot = clonal_solutions %>% 
    ggplot(aes(x = k1, y = loglikelihood)) +
    ggplot2::geom_col(fill = '#2274A5')+ my_ggplot_theme() +#+ scale_fill_manual(c('best'='#E574BC', ''= '#2274A5')) + my_ggplot_theme() 
    labs(x= 'CCF') +
    ggplot2::ggtitle('Clonal solutions')
}

#' Plot summarising the results of patch for a single segment (best solution, all subclonal solutions, all clonal solutions)
#'
#' @param x a cnaqc object
#' @param seg_id the segment id of the segment of interest
#'
#' @return
#' @export
#'
#' @examples
patch_plot = function(x,seg_id){
  ggarrange(plotlist= list(plot_patch_best_solution(x, seg_id), plot_patch_all_solutions_subclonal(x,seg_id),plot_patch_all_solutions_clonal(x,seg_id) ), ncol=3)
}