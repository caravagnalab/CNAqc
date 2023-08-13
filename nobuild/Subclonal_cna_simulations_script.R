devtools::load_all('~/Documents/Documents/GitHub/CNAqc/R')

options = expand.grid(purity = seq(.1, 1, .1), coverage = seq(30, 100, 10), ccf = seq(.1, .9, .1), replicates = seq(1, 10, 1))

overall_res = lapply(seq(1, 2,#length(options$purity),
                         by =1), function(i){
  x = simulate_sample(t_subclone_mrca = 5, 
                      cna_rate = 5,
                      mu= 6e-8, w= 10, 
                      coverage = options[i,]$coverage, purity= options[i,]$purity,
                      sub_cnas_rate= 7, t_f= 10,
                      ccf = options[i,]$ccf, bin_size= 50000)
  
  mutations = x$mutations %>% mutate(Major= as.integer(Major),  minor = as.integer(minor)) %>% select(!segment_id)#%>% mutate(ref = 'A', alt = 'T')
  cna = x$cna %>% mutate(Major= as.integer(Major),  minor = as.integer(minor)) %>% mutate(CCF=1) %>% select(!segment_id)
  snps = x$snps %>% select(!segment_id)
  ref_genome = 'hg19'
  purity = x$purity
  x = init(mutations=mutations, cna= cna, snps = snps, purity= purity, sample = "simulation", ref = ref_genome)
  
  patched_x = patch(x, all_solutions=TRUE)
  dir = paste0('Coverage_',options[i,]$coverage, '_Purity_', options[i,]$purity*10, '_CCF_',ccf = options[i,]$ccf*10)
  
  if (!file.exists(paste0('../nobuild/simulations/', dir))){
    dir.create(paste0('../nobuild/simulations/', dir))
  }
  
  original_res = x$cna %>% mutate('k1'= paste0(Major, ':', minor), k2 = ifelse(Major_2==0, '', paste0(Major_2, ':', minor_2))) %>%
    select(k1, k2, CCF, segment_id)
  
  new_res = patched_x$patch_best_solution %>% 
    mutate('k1_best'= k1, 'k2_best' = k2, 'CCF_best'= ccf_1) %>%
    select(k1_best, k2_best, CCF_best, segment_id)
  
  merged_res = merge(original_res, new_res) %>% 
    mutate(k_res = (k1==k1_best & k2==k2_best), ccf_res = (CCF==CCF_best), all = (k1==k1_best & k2==k2_best & CCF==CCF_best) )
  res_plot= merged_res %>% select(segment_id, k_res, ccf_res, all) %>% as_tibble() %>% reshape2::melt(id = 'segment_id') %>% ggplot() +
    geom_bar(aes(x=variable, fill=value)) + theme_bw() + labs(x='', y='')+scale_fill_manual(values = c('TRUE'='forestgreen', 'FALSE'='indianred'))+
    scale_x_discrete(labels=c('matched karyotype', 'matched ccf', 'match'))+
    theme(
      axis.text.x = element_text(angle = 45, hjust=1)
    )
  
  segment_ids = patched_x$patch_best_solution %>% pull(segment_id) %>% unique()
  
  for (j in 1:length(segment_ids)){
    #print(segment_ids[i])
    pp = patch_plot(patched_x, segment_ids[j])
    ggsave(pp, filename= paste0('../nobuild/simulations/', dir,'/', segment_ids[j], '.png'), width = 20, height = 10)
  }
  
  
  patched_file_name = paste0('../nobuild/simulations/', dir, '/', options[i,]$replicates, '.rds')
  csv_file_name = paste0('../nobuild/simulations/', dir, '/', options[i,]$replicates, '.csv')
  
  saveRDS(patched_x, patched_file_name)
  write.csv(merge(original_res, new_res), file= csv_file_name)
  ggsave(res_plot, filename=paste0('../nobuild/simulations/', dir, '/res_plot_', options[i,]$replicates, '.png'))
  
  summary_table = merge(original_res, new_res) %>% mutate(purity = options[i,]$purity, 
                                                          coverage = options[i,]$coverage, 
                                                          replicate = options[i,]$replicate,
                                                          ccf = options[i,]$ccf)
  return(summary_table)
})

overall_res = Reduce(rbind, overall_res) 
write.csv(overall_res, file= '../nobuild/simulations/simulation_results.csv')

overall_res =overall_res %>%
  mutate(k_res = (k1==k1_best & k2==k2_best), ccf_res = (CCF==CCF_best), all = (k1==k1_best & k2==k2_best & CCF==CCF_best),
         sim_id = paste0('Purity_', purity,'_Coverage_',coverage,'_CCF_', ccf))

res_plot= overall_res %>% select(sim_id, k_res, ccf_res, all) %>% as_tibble() %>% 
  reshape2::melt(id = 'sim_id') %>% ggplot() +
  geom_bar(aes(x=variable, fill=value)) + theme_bw() + labs(x='', y='')+
  scale_fill_manual(values = c('TRUE'='forestgreen', 'FALSE'='indianred'))+
  scale_x_discrete(labels=c('matched karyotype', 'matched ccf', 'match'))+
  theme(
    axis.text.x = element_text(angle = 45, hjust=1)
  )+ facet_wrap(~ sim_id, scales= 'free') 

ggsave(res_plot, filename=paste0('../nobuild/simulations/res_plot.png'), height = 40, width = 40)

# Single Simulation
# 
# x = simulate_sample(t_subclone_mrca = 5, cna_rate = 5,
#                     mu= 6e-8, w= 10, coverage = 100, purity=1,
#                     sub_cnas_rate= 7, t_f= 10, ccf = .5, bin_size= 50000)
# # Check simulation
# plot_snps(x, s= 10)
# plot_snps(x, s= 10, what = 'DR')
# x$mutations %>% ggplot() + geom_histogram(aes(x= VAF))+ facet_wrap(~segment_id) + theme_bw() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
# relative_to_absolute_coordinates(x, x$cna) %>% mutate(Major= as.integer(Major),  minor = as.integer(minor)) %>% ggplot() + 
#   geom_segment(aes(x= from, xend=to, y= Major+.1, yend= Major+.1), color = 'indianred')+
#   geom_segment(aes(x= from, xend=to, y= minor -.1, yend= minor -.1), color = 'steelblue')+
#   theme_bw()+ labs(x = 'absolute genome coordinate', y = 'allele multiplicity (clonal)')
# 
# mutations = x$mutations %>% mutate(Major= as.integer(Major),  minor = as.integer(minor)) %>% select(!segment_id)#%>% mutate(ref = 'A', alt = 'T')
# cna = x$cna %>% mutate(Major= as.integer(Major),  minor = as.integer(minor)) %>% mutate(CCF=1) %>% select(!segment_id)
# snps = x$snps %>% select(!segment_id)
# ref_genome = 'hg19'
# purity = x$purity
# 
# new_cnaqc_obj = init(mutations=mutations, cna= cna, snps = snps, purity= purity, sample = "simulation", ref = ref_genome)
# 
# new_cnaqc_obj = patch(new_cnaqc_obj, all_solutions=TRUE)
# saveRDS(new_cnaqc_obj, '../nobuild/simulated_sample.rds')
# 
# segment_ids = new_cnaqc_obj$patch_best_solution %>% pull(segment_id) %>% unique()
# 
# for (i in 1:length(segment_ids)){
#   print(segment_ids[i])
#   pp = patch_plot(new_cnaqc_obj, segment_ids[i])
#   ggsave(pp, filename= paste0('../nobuild/simulation_trial/', segment_ids[i], '.png'), width = 20, height = 10)
# }
# 
# original_res = x$cna %>% mutate('k1'= paste0(Major, ':', minor), k2 = ifelse(Major_2==0, '', paste0(Major_2, ':', minor_2))) %>%
#   select(k1, k2, CCF, segment_id)
# 
# new_res = new_cnaqc_obj$patch_best_solution %>% 
#   mutate('k1_best'= k1, 'k2_best' = k2, 'CCF_best'= ccf_1) %>%
#   select(k1_best, k2_best, CCF_best, segment_id)
# 
# merged_res = merge(original_res, new_res) %>% 
#   mutate(k_res = (k1==k1_best & k2==k2_best), ccf_res = (CCF==CCF_best), all = (k1==k1_best & k2==k2_best & CCF==CCF_best) )
# merged_res %>% select(segment_id, k_res, ccf_res, all) %>% as_tibble() %>% reshape2::melt(id = 'segment_id') %>% ggplot() +
#   geom_bar(aes(x=variable, fill=value)) + theme_bw() + labs(x='', y='')+scale_fill_manual(values = c('TRUE'='forestgreen', 'FALSE'='indianred'))+
#   scale_x_discrete(labels=c('matched karyotype', 'matched ccf', 'match'))+
#   theme(
#     axis.text.x = element_text(angle = 45, hjust=1)
#   )
# 
# 
