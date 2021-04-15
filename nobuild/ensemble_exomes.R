library(biomaRt)

ensembl = useEnsembl(
  biomart="ensembl",
  dataset="hsapiens_gene_ensembl",
  # host = "https://grch37.ensembl.org" hg19 (remove to get GRCh38)
  )

chr_ex_gen_all = NULL

for (chr in c(1:22, 'X', 'Y'))
{
  cli::cli_h1(chr)

  exones_chr <- getBM(
    attributes = c(
      'ensembl_exon_id',
      'chromosome_name',
      'exon_chrom_start',
      'exon_chrom_end'
    )
    ,
    filters =
      'chromosome_name',
    values = chr,
    mart = ensembl
  )
  exones_chr =  exones_chr  %>% as_tibble()

  cli::cli_h3(chr)

  gene_exones <- getBM(
    attributes = c('hgnc_symbol',
                   'ensembl_exon_id')
    ,
    filters =
      'chromosome_name',
    values = chr,
    mart = ensembl
  )
  gene_exones = gene_exones %>% as_tibble() %>% filter(hgnc_symbol != '')

  chr_ex_gen = gene_exones %>%
    left_join(exones_chr) %>%
    distinct() %>%
    mutate(chromosome_name = paste(chromosome_name))

  cli::cli_alert(object.size(chr_ex_gen))

  chr_ex_gen %>% print

  saveRDS(chr_ex_gen, file = paste0(chr, '.rds'))
  chr_ex_gen_all = chr_ex_gen_all %>% bind_rows(chr_ex_gen)
}

chr_ex_gen_all = chr_ex_gen_all[complete.cases(chr_ex_gen_all),]
chr_ex_gen_all$chromosome_name = paste0('chr', chr_ex_gen_all$chromosome_name)
chr_ex_gen_all = chr_ex_gen_all %>%
  rename(
    gene = hgnc_symbol,
    chr = chromosome_name,
    from = exon_chrom_start,
    to = exon_chrom_end
  ) %>%
  dplyr::select(gene, chr, from, to, ensembl_exon_id)

# exones_GRCh38 = exones_GRCh38 %>% distinct(gene, chr, from, to)
exones_hg19 = chr_ex_gen_all %>% distinct(gene, chr, from, to)

# usethis::use_data(exones_GRCh38, overwrite = TRUE)
usethis::use_data(exones_hg19, overwrite = TRUE)

annotate_drivers()

