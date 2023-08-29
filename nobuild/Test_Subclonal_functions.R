devtools::load_all('~/Documents/Documents/GitHub/CNAqc/R')

x = readRDS('../data/Test_subclonalCNA_organoid.rds')

mutations = x$mutations
cna = x$cna
ref_genome = x$reference_genome
purity = x$purity

## SNPs dataframe
file_dir = '/Users/aliceantonello/Library/CloudStorage/Dropbox/Organoids_Accelerator/data/Sequenza/PDO_74'
r= load(paste0(file_dir,'/74_PDO_2_sequenza_extract.RData'))
r = get(r)  

BAF_DR = extract_sequenza_baf_dr(r)
snps = BAF_DR$binned

new_cnaqc_obj = init(mutations=mutations, cna= cna, snps = snps, purity= purity, sample = "PDO74", ref = ref_genome)
#print(new_cnaqc_obj)

#inspect_segment(new_cnaqc_obj)
#baf_plot= plot_snps(new_cnaqc_obj)
#dr_plot= plot_snps(new_cnaqc_obj, what = 'DR')
#my_baf_dr_plot = baf_plot / dr_plot

new_cnaqc_obj = patch(new_cnaqc_obj, all_solutions=TRUE, baf_coef = 15, dr_coef = 10, vaf_coef= 20)
saveRDS(new_cnaqc_obj, '../nobuild/new_cnaqc_object_organoid.rds')

segment_ids = new_cnaqc_obj$patch_best_solution %>% 
  mutate(from = stringr::str_split(segment_id, pattern= ':', simplify = T)[, 2], 
         to = stringr::str_split(segment_id, pattern= ':', simplify = T)[, 3]) %>% 
  mutate(from = as.integer(from), to= as.integer(to)) %>%
  mutate(l= to-from) %>% filter(l>3e7) %>% pull(segment_id) %>% unique()

for (i in 1:length(segment_ids)){
  print(segment_ids[i])
  pp = patch_plot(new_cnaqc_obj, segment_ids[i])
  ggsave(pp, filename= paste0('../nobuild/PDO74_results/plots1/', segment_ids[i], '.png'), width = 20, height = 10)
}


new_vs_old_res = function(x, selected_segment){
  
  seg1_muts= x$mutations %>% filter(segment_id==selected_segment)
  seg1_snps = x$snps %>% filter(segment_id==selected_segment)
  
  my_new_plot_baf_exp_vs_obs = function(SNPs, k1 = '1:1', k2 = '', purity = 1, model = '', sim_color='steelblue', ccf=0, s=5){
    
    n_points = length(SNPs$chr)
    n_snps = SNPs %>% pull(N.BAF) %>% sum()
    
    if (k2 == ''){
      
      expected_baf_d = expected_baf_ditribution(n= n_snps, k1= k1, k2='', 
                                                purity =  purity, ccf=1, g1='', g2='', n_points)
    }else{
      
      discrete_ccf = seq(from=0, to= 1, by= .05)
      temp_pur= discrete_ccf[which(abs(discrete_ccf - purity) == min(abs(discrete_ccf - purity)))]
      
      genotypes = get_models(karyotype_combination = paste0(k1,'-' ,k2),
                             purity = temp_pur, ccf = ccf) %>% filter(model_id == model) %>% 
        select(genotype_1,genotype_2) 
      g1 = genotypes %>% filter(!is.na(genotype_1)) %>% pull(genotype_1) %>% unique()
      g2 = genotypes %>% filter(!is.na(genotype_2)) %>% pull(genotype_2) %>% unique()
      expected_baf_d = expected_baf_ditribution(n= n_snps, k1= k1, k2= k2, 
                                                purity =  purity, ccf= ccf, g1, g2, n_points)
    }
    
    e_baf_d = SNPs %>% 
      select(from, to) %>% mutate(exp_baf = expected_baf_d) %>% as_tibble() %>% mutate(from=as.integer(from), to=as.integer(to))
    
    # the observed baf distribution is constant 
    observed_baf_d = SNPs %>% select(from, to, BAF)
    
    ggplot() +
      geom_segment(data= e_baf_d, aes(x = from, xend= to, y= exp_baf, yend = exp_baf), color= sim_color, size = s)+
      #geom_density_2d(data= e_baf_d, aes(x = from, y= exp_baf))+
      geom_segment(data = observed_baf_d %>% as_tibble(), aes(x = from, xend= to, y= BAF, yend = BAF), color= 'black', size = s)+
      theme_bw() + labs(x= '', y= 'BAF') + ylim(0,1)
  }
  my_new_plot_dr_exp_vs_obs = function(SNPs, k1 = '1:1', k2 = '', purity = 1, model = '', sim_color='steelblue', ccf=0, s=5, ploidy=2){
    
    n_points = length(SNPs$chr)
    n_snps = SNPs %>% pull(N.BAF) %>% sum()
    
    if (k2==''){
      expected_dr_d = expected_dr_distribution(n= n_snps, k1= k1, k2='', 
                                               purity =  purity, ccf=1, ploidy = ploidy, n_points)
    }else{
      expected_dr_d = expected_dr_distribution(n= n_snps, k1= k1, k2=k2, 
                                               purity =  purity, ccf=ccf, 
                                               ploidy = ploidy, n_points)
    }
    
    e_dr_d = SNPs %>% 
      select(from, to) %>% mutate(exp_dr = expected_dr_d) %>% as_tibble() %>% mutate(from=as.integer(from), to=as.integer(to))
    # the observed baf distribution is constant 
    observed_dr_d = SNPs %>% select(from, to, DR)
    
    ggplot() +
      geom_segment(data= e_dr_d, aes(x = from, xend= to, y= exp_dr, yend = exp_dr), color= sim_color, size=s)+
      geom_segment(data = observed_dr_d %>% as_tibble(), aes(x = from, xend= to, y= DR, yend = DR), color= 'black',size=s)+
      theme_bw() + labs(x= '', y= 'DR')
  }
  my_new_plot_vaf_exp_vs_obs = function(SNVs, k1 = '1:1', k2 = '', purity = 1, model = '', sim_color='steelblue', ccf=0){
    
    n_points = length(SNVs$chr)
    coverage = mean(SNVs$DP)
    
    if (k2==''){
      peaks = get_clonal_peaks(k1, purity)
    }else{
      
      discrete_ccf = seq(from=0, to= 1, by= .05)
      temp_pur= discrete_ccf[which(abs(discrete_ccf - purity) == min(abs(discrete_ccf - purity)))]
      
      peaks = get_models(karyotype_combination = paste0(k1,'-' ,k2),
                         purity = temp_pur, ccf = ccf) %>% filter(model_id == model) %>% pull(peak)
    }
      
    vd = vaf_distribution(n_points, coverage, peaks)
    ggplot(vd) +
      geom_histogram(aes(x= nv/dp), fill= sim_color, alpha=.5)+
      geom_histogram(data = SNVs, aes(x= VAF), fill= 'black', alpha=.5)+
      theme_bw()+labs(y='', x='VAF') + xlim(0,1)
  }
  
  major_old = x$cna %>% filter(segment_id== selected_segment) %>% pull(Major)
  minor_old = x$cna %>% filter(segment_id== selected_segment) %>% pull(minor)
  k1_old = paste0(major_old, ':', minor_old)
  
  old_baf = my_new_plot_baf_exp_vs_obs(seg1_snps, purity= x$purity,
                             k1 = k1_old)
  old_dr = my_new_plot_dr_exp_vs_obs(seg1_snps, purity= x$purity, ploidy= x$ploidy,
                            k1 = k1_old)
  old_vaf = my_new_plot_vaf_exp_vs_obs(seg1_muts, purity= x$purity,
                             k1 = k1_old)
  
  old_title = paste0('Clonal segment with karyotype ', k1_old)
  
  st = 'AAACC
        BBBCC'
  
  old_res_plot = patchwork::wrap_plots(old_baf, old_dr, old_vaf, design = st) + patchwork::plot_annotation(title = old_title)
  
  k1_best = x$patch_best_solution[1,]$k1
  k2_best = x$patch_best_solution[1,]$k2
  model_best = x$patch_best_solution[1,]$model
  ccf_best = x$patch_best_solution[1,]$ccf
  
  if (ccf_best==1){
    new_title = paste0('Clonal segment with karyotype ', k1_best)
  }else{
    new_title = paste0('Sublonal segment with karyotypes ', k1_best, ' , ', k2_best, ' ccf ', ccf_best)
  }
  
  
  new_baf = my_new_plot_baf_exp_vs_obs(seg1_snps, purity= x$purity, k1 = k1_best, 
                             k2 = k2_best,  
                             model = model_best, 
                             sim_color='darkorange', 
                             ccf= ccf_best)
  new_dr = my_new_plot_dr_exp_vs_obs(seg1_snps, purity= x$purity, 
                            ploidy= x$ploidy, k1 = k1_best, k2 = k2_best,  model = model_best, 
                            sim_color='darkorange', ccf=ccf_best)
  new_vaf = my_new_plot_vaf_exp_vs_obs(seg1_muts, purity= x$purity, k1 = k1_best, 
                             k2 = k2_best,  model = model_best, sim_color='darkorange', ccf=ccf_best)
  
  new_res_plot = patchwork::wrap_plots(new_baf, new_dr, new_vaf, design = st) + patchwork::plot_annotation(title = new_title)
  
  final_plot = ggarrange(old_res_plot, new_res_plot, ncol=2)
  final_plot
}

segments = new_cnaqc_obj$patch_best_solution %>% pull(segment_id) %>% unique()

lapply(segments, function(i){
  res_plot= new_vs_old_res(new_cnaqc_obj,i)
  ggsave(res_plot, filename= paste0('../nobuild/PDO74_results/plots2/', i, '.png'), width = 20, height = 10)
})


#seg2 = new_cnaqc_obj$mutations %>% filter(segment_id=="chr1:624505:124997413:1:1:1")





