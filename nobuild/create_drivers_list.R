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
  



