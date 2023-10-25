# Create list of all processes ####
karyotypes_1 <- c("1:0", "1:1", "2:0", "2:1", "2:2")
karyotypes_2 <- c("1:0", "1:1", "2:0", "2:1", "2:2")

k1_used <- c()
total_list_of_processes <- list()
for (k1 in karyotypes_1) {
  print(paste0('karyotype 1 = ',k1))
  for (k2 in karyotypes_2) {
    if (sum(karyotype_distances(k1, k2)) <= 2 && !(k2 %in% k1_used)) {
      if (sum(karyotype_distances(k1, k2)) <= 2) {
        print(paste0('     karyotype 2 = ',k2))
        list_of_processes <- get_evolutionary_processes('1:1', k1, k2)
        model_id <- paste(k1, k2, sep='-')
        total_list_of_processes[[model_id]] <- list_of_processes
      }
    }
  }
  k1_used <- c(k1_used, k1)
}
#saveRDS(total_list_of_processes, './segment_scripts//NEW_total_list_of_processes.rds')
total_list_of_processes


#list_of_processes <- readRDS("./segment_scripts/extra/NEW_total_list_of_processes.rds")
clonal_karyotypes <- c('1:0', '1:1', '2:0', '2:1', '2:2')

obtain_models <- function(total_list_of_processes){
  a <- list()
  names <- c()
  
  for (n in names(total_list_of_processes)) {
    
    names <- c(names, n)
    
    c <- list()
    kr <- c()
    
    processes <- total_list_of_processes[[n]]
    karyotypes <- strsplit(n[1], "-")[[1]]
    k1 <- karyotypes[1]
    k2 <- karyotypes[2]
    
    for (ccf in seq(.1, .9, by = 0.1)){
      for (pr in seq(0, 1, by = 0.01)){
        
        lab <- paste(ccf,pr,sep='-')
        print(lab)
        
        dd <- expectations_subclonal(list_of_processes = processes, CCF_1 = ccf, purity = pr)
        d <- dd %>% select(peak, model_id, genotype_1, genotype_2)
        
        c <- append(c, list(d))
        kr <- c(kr, lab)
      }
    }
    
    names(c) <- kr
    a <- append(a, list(c))
  }
  names(a) <- names
  return(a)
}
res <- obtain_models(total_list_of_processes)
saveRDS(res, '/data/model_ccf_purity.rds')
#mm <- readRDS('./segment_scripts/model_ccf_purity.rds')


