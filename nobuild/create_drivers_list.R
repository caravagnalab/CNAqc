DriverDBv3 = readr::read_tsv("./mutation_download_tab.txt")

ntools = unique(DriverDBv3$tool)

DriverDBv3 = DriverDBv3 %>%
  apply(MARGIN = 1, FUN = function(w) {
    data.frame(
      gene = strsplit(w['driver_gene'], ', '),
      cancer = w['cancer_project'],
      tool = w['tool'],
      stringsAsFactors = FALSE
    )
  }) %>%
  Reduce(f = bind_rows) %>%
  as_tibble() %>%
  arrange(driver_gene) %>%
  filter(driver_gene != "")

DriverDBv3 = DriverDBv3 %>%
  group_by(driver_gene, cancer) %>%
  summarise(n_tools = n()) %>%
  arrange(desc(n_tools)) %>%
  filter(n_tools > 3) %>%
  arrange(driver_gene)

# DriverDBv3 %>%
#   filter(driver_gene %in% DriverDBv3_to_keep$driver_gene) %>%
#   group_by(driver_gene, cancer) %>%
#   summarise(tools = paste(tool, collapse ="|"))


usethis::use_data(DriverDBv3, overwrite = TRUE)

# exones_GRCh38 =  readr::read_tsv("~/GR_txt.tsv", col_names = F)
# colnames(exones_GRCh38) = c("chr", "from", "to")
#
# gene_coordinates_GRCh38
#
# i = 100
# what_exones_chr = exones_GRCh38 %>%
#   filter(chr == gene_coordinates_GRCh38$chr[i]) %>%
#   filter(from >= gene_coordinates_GRCh38$from[i]) %>%
#   filter(to <= gene_coordinates_GRCh38$to[i])
#
# (exones_GRCh38$from >= gene_coordinates_GRCh38$from[i]) &
#   (exones_GRCh38$to <= gene_coordinates_GRCh38$to[i])
#
# exones_GRCh38[what, ]



