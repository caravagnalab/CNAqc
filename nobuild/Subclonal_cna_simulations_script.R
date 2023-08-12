devtools::load_all('~/Documents/Documents/GitHub/CNAqc/R')

x = simulate_sample(t_subclone_mrca = 5, cna_rate = 5,
                    mu= 6e-8, w= 10, coverage = 100, purity=1,
                    sub_cnas_rate= 5, t_f= 10, ccf = .5, bin_size= 50000)

x$snps %>% ggplot() + geom_point(aes(x= from, y= BAF)) + facet_wrap(~segment_id, scales='free')
x$snps %>% ggplot() + geom_point(aes(x= from, y= DR)) + facet_wrap(~segment_id, scales='free')
x$mutations %>% ggplot() + geom_histogram(aes(x= VAF))+ facet_wrap(~segment_id)
relative_to_absolute_coordinates(x, x$cna) %>% ggplot() %>% geom_segment(aes(x= from, xend=to, y= Major, yend= Major, color= model_id))

mutations = x$mutations %>% mutate(Major= as.integer(Major),  minor = as.integer(minor)) #%>% mutate(ref = 'A', alt = 'T')
cna = x$cna %>% mutate(Major= as.integer(Major),  minor = as.integer(minor)) # %>% mutate(CCF=1)
snps = x$snps 
ref_genome = 'GRCh37'
purity = x$purity

new_cnaqc_obj = init(mutations=mutations, cna= cna, snps = snps, purity= purity, sample = "simulation", ref = ref_genome)
#print(new_cnaqc_obj)

#plot_snps(new_cnaqc_obj, what = 'BAF')
#plot_snps(new_cnaqc_obj, what = 'DR')

new_cnaqc_obj = patch(new_cnaqc_obj, all_solutions=TRUE, preselect = TRUE)
#saveRDS(new_cnaqc_obj, 'new_cnaqc_object_organoid.rds')

segment_ids = new_cnaqc_obj$patch_best_solution %>% pull(segment_id)


for (i in 1:length(segment_ids)){
  pp= patch_plot(new_cnaqc_obj, segment_ids[i])
  ggsave(pp, filename= paste0('../nobuild/simulation_trial/', segment_ids[i], '.png'))
}




