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
#extract_sequenza_baf_dr(r)
#x$BAF_DR = r

BAF_DR = extract_sequenza_baf_dr(r)
snps = BAF_DR$binned

new_cnaqc_obj = init(mutations=mutations, cna= cna, snps = snps, purity= purity, sample = "PDO74", ref = ref_genome)
print(new_cnaqc_obj)


new_cnaqc_obj = patch(new_cnaqc_obj, all_solutions=TRUE, preselect = TRUE)
saveRDS(new_cnaqc_obj, '../nobuild/new_cnaqc_object_organoid.rds')

segment_ids = new_cnaqc_obj$patch_best_solution %>% pull(segment_id) %>% unique()
for (i in 1:length(segment_ids)){
  print(segment_ids[i])
  pp = patch_plot(new_cnaqc_obj, segment_ids[i])
  ggsave(pp, filename= paste0('../nobuild/PDO74_results/', segment_ids[i], '.png'), width = 20, height = 10)
}



