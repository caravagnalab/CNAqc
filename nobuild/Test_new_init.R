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

plot_snps(new_cnaqc_obj, what = 'BAF')
plot_snps(new_cnaqc_obj, what = 'DR')


