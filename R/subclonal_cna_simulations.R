# Simulate a CNAqc object

# SNVs (clonal / tail / subclonal)
# SNPs 
# CNAs 
# purity
# ploidy : weighted average of karyotype

# Times : t0 = 0, t_subclonal_mrca = 5, t_final = 10

# 1. SIMULATE THE CLONAL CNAs 
# 1.1 Fix the amount of time from 0 to appereance of first subclone (drawn from a uniform [0, t_max])
# 1.2 Fix the rate of CNA (cna_rate)
# 1.3 Draw a number of CNAs from a Poisson distribution with rate (cna_rate) and 
#   the corresponding time of occurrence from a uniform
# 1.4. Draw their length and position from a uniform distribution:
#         - chromosome -> length -> from/to -> model (final karyotypes and alleles composition)

sample_karyo = function(){
  k = c('1:0', '2:1', '2:0', '2:2')
  k[sample.int(4, 1)]
}
sample_chr = function(filtered_chr=NULL){
  if (is.null(filtered_chr)){chromosomes = c(paste0('chr', 1:22), 'chrX', 'chrY')}
  else{chromosomes = filtered_chr}
  chromosomes[sample.int(length(chromosomes), 1)]
}
sample_allele = function(g){
  have_allele = FALSE
  available_alleles = c('A1', 'A2', 'B1', 'B2')
  while (!have_allele){
    a = available_alleles[sample.int(4, 1)]
    if (grepl(a, g, fixed = TRUE)){
      have_allele=TRUE
    }
  }
  a
}
get_chr_length = function(c){
  CNAqc:::get_reference("GRCh37") %>% filter(chr == c) %>% pull(length)
}
sample_from_to = function(c){
  chr_l = get_chr_length(c)
  l= runif(1, 1e7, chr_l)
  from = sample.int(chr_l-l, 1)
  to = from + l
  
  return(c(from, to))
}
sample_genotype = function(k){
  if (k=='1:0'){
    g = c('A1', 'B1')
    return(g[sample.int(length(g), 1)])
    }
  if (k=='2:1'){
    g = c('A1A2B1', 'A1B1B2')
    return(g[sample.int(length(g), 1)])
  }
  if (k=='2:0'){
    g = c('A1A2', 'B1B2')
    return(g[sample.int(length(g), 1)])
  }
  if (k=='2:2'){
    g = c('A1A2B1B2')
    return(g)
  }
  
}
sample_model = function(k, purity, ccf){
  kc = paste0('1:1-',k)
  if (kc %in% names(get_models())){
    
    subclonal_cn_df = get_models(karyotype_combination= kc, purity= purity, ccf= ccf) 
    
    if (length(subclonal_cn_df$model_id) >=1){
      model_ids = subclonal_cn_df %>% pull(model_id) %>% unique()
      target_model = model_ids[sample.int(length(model_ids), 1)]
      
      return(subclonal_cn_df %>% filter(model_id == target_model))
    }else{
      return('invalid karyotype')
    }
  }else{
    return('invalid karyotype')
  }
  
}
fill_the_gaps = function(simulated_cnas){
  
  chromosomes = c(paste0('chr', 1:22), 'chrX', 'chrY')
  filled_cna = lapply(chromosomes, function(c){
    l = get_chr_length(c)
    chromosomes_c = simulated_cnas %>% pull(chr)
    if (!(c %in% chromosomes_c)){
      new_row = data.frame('chr' = c, 'from'=0, 'to'= l, 'karyotype'= '1:1', 'genotype'= 'A1B1', 'times'= 0)
    }else{
      simulated_cnas_c = simulated_cnas %>% filter(chr == c)
      n_cna_c = length(simulated_cnas_c %>% pull(chr))
      new_row = data.frame()
      for (i in 1:n_cna_c){
        if (i==1){
          new_row_11 = data.frame('chr' = c, 'from'=0, 'to'= simulated_cnas_c[1,]$from, 'karyotype'= '1:1', 'genotype'= 'A1B1', 'times'= 0)
          new_row_12 = data.frame('chr' = c, 'from'=simulated_cnas_c[1,]$to, 'to'= l, 'karyotype'= '1:1', 'genotype'= 'A1B1', 'times'= 0)
          new_row_i = rbind(new_row_11, new_row_12)
        }else{
          if (i == n_cna_c){
            new_row_i = data.frame('chr' = c, 'from'= simulated_cnas_c[n_cna_c,]$to, 'to'= l, 'karyotype'= '1:1', 'genotype'= 'A1B1', 'times'= 0)
          }else{
            new_row_i = data.frame('chr' = c, 'from'= simulated_cnas_c[i-1,]$to, 'to'= simulated_cnas_c[i,]$from, 'karyotype'= '1:1', 'genotype'= 'A1B1', 'times'= 0)
          }
        }
        new_row = rbind(new_row_i, new_row)
        
      }
    }
    new_row
  })
  
  filled_cna = Reduce(rbind, filled_cna) 
  filled_cna = rbind(simulated_cnas, filled_cna)
  filled_cna %>% arrange(chr)
}
has_subclonal_event = function(cna, chrom){
  (length(cna %>% filter(chr==chrom) %>% pull(chr)) > 1)
}
my_pareto_betabin <- function(N){
  tail.cutoff <- .5
  shape <- runif(1,1,3) # shape of the pareto depends on mutation rate
  scale <- .05
  tails <- sads::rpareto(N, shape = shape, scale = scale)
  tails <- tails[tails <= tail.cutoff]
  return(tails) # return vaf of subclonal mutations
}
compute_approx_ploidy = function(cnas){
  ploidy = cnas %>% filter(CCF==1) %>% mutate(total_cn = as.integer(Major) + as.integer(minor)) %>% pull(total_cn) %>% mean()
  ploidy
}

simulate_clonal_cnas = function(t_subclone_mrca = 5, cna_rate = 5){
  n_cnas = rpois(1, cna_rate)
  times = runif(n_cnas, 0, t_subclone_mrca) %>% sort()
  
  clonal_cnas = data.frame()
  for (i in 1:n_cnas){
    karyotypes = sample_karyo()
    chrs = sample_chr()
    form_to = sample_from_to(chrs)
    from = form_to[1]
    to = form_to[2]
    genotype = sample_genotype(karyotypes) 
    cc = data.frame('chr' = chrs, 'from'=from, 'to'= to, 'karyotype'= karyotypes, 'genotype'= genotype)
    clonal_cnas = rbind(cc, clonal_cnas)
  }
  clonal_cnas = clonal_cnas %>% mutate('times'= times) %>% fill_the_gaps() %>% rowwise() %>% 
    mutate('Major'= strsplit(karyotype, ':')[[1]][1], 'minor' = strsplit(karyotype, ':')[[1]][2]) %>%
    mutate(segment_id = paste0(chr,':' ,from,':' ,to, ':',Major, ':',minor, ':1'))
}

# 2. SIMULATE CLONAL SNVs 
# 2.1 Place the clonal snvs on each segment :
#         - the older ones with higher multiplicities will accumulate with rate depending on 
#           mutation rate and growth rate, over a period of time going from 0 to time of segment formation
#         - the newer ones with rate depending on a period going from the time of segment formation to subclone birth

generate_clonal_snvs = function(seg_id, karyotype, chromosome, al='A1B1',
                                mu= 1e-8, w= 2, t_subclone_mrca = 5, 
                                coverage = 100, purity=1, t_cna){
  
  if (karyotype == '1:0'){
    n_snvs = rpois(1, mu*w*t_subclone_mrca*get_chr_length(chromosome)) 
  }
  if (karyotype == '1:1'){
    n_snvs = rpois(1, 2*mu*w*t_subclone_mrca*get_chr_length(chromosome)) 
  }
  if (karyotype == '2:0'){
    n_snvs = rpois(1, mu*w*t_cna*get_chr_length(chromosome)) + rpois(1, 2*mu*w*(t_subclone_mrca-t_cna)*get_chr_length(chromosome))
  }
  if (karyotype == '2:1'){
    n_snvs = rpois(1, 2*mu*w*t_cna*get_chr_length(chromosome)) + rpois(1, 3*mu*w*(t_subclone_mrca-t_cna)*get_chr_length(chromosome))
  }
  if (karyotype == '2:2'){
    n_snvs = rpois(1, 2*mu*w*t_cna*get_chr_length(chromosome)) + rpois(1, 4*mu*w*(t_subclone_mrca-t_cna)*get_chr_length(chromosome))
  }
  
  dp = rpois(n_snvs, coverage)
  peak_mean = get_clonal_peaks(karyotype, purity)
  p = rbeta(n_snvs, round(20*peak_mean), round(20*(1-peak_mean)))
  nv = rbinom(n_snvs, dp, p)
  from = sample.int(get_chr_length(chromosome)-1, n_snvs)
  to = from +1
  allele = c()
  for (i in 1:n_snvs){allele[i]=sample_allele(al)}
  #allele = rep(sample_allele('A1B1'), n_snvs)
  maj = strsplit(karyotype, ':')[[1]][1]
  minor = strsplit(karyotype, ':')[[1]][2]
  snvs_df = data.frame('chr' = chromosome, 'from'=from, 'to'= to, 
                       'Major'= maj,'minor'=minor,'NV'=nv,'DP'=dp, 'VAF'=nv/dp,
                       'allele'= allele, 'segment_id'= seg_id)
  
}

simulate_clonal_snvs_df = function(cnas, mu= 1e-8, w= 2, t_subclone_mrca = 5, coverage = 100, purity=1){
  segments = cnas$segment_id
  clonal_snvs_df = lapply(segments, function(x){
    chromosome = strsplit(x, ':')[[1]][1]
    
    snvs_df = data.frame()
    cna_s = cnas %>% filter(segment_id ==x)
    
    snvs_df_i = generate_clonal_snvs(seg_id=x, karyotype= cna_s$karyotype,
                                     chromosome=chromosome, al= cna_s$genotype,
                                     mu=mu, w= w, t_subclone_mrca = t_subclone_mrca,
                                     coverage = coverage, purity=purity, t_cna=cna_s$times)
    snvs_df = rbind(snvs_df_i, snvs_df)
  })
  clonal_snvs_df = Reduce(rbind, clonal_snvs_df)
  clonal_snvs_df
  
}


# 3. SIMULATE THE SUBCLONAL CNAs
# 3.1 Fix the final CCF
# 3.2 Draw a number of subclonal CNAs from (lam_subclonal)
# 3.3 Draw their length and position from a uniform distribution : chromosome -> length -> from/to 
# 3.4 Based on which segment they fall on, determine the model (allele composition)

simulate_subclonal_cnas = function(cnas, sub_cnas_rate= 3, t_subclone_mrca = 5, t_final= 10, ccf = .5, purity =1){
  
  n_sub_cnas = rpois(1, sub_cnas_rate)
  chromosomes = c(paste0('chr', 1:22), 'chrX', 'chrY')
  chromosomes_without_cnas = chromosomes[!sapply(chromosomes, has_subclonal_event, cna=cnas)]
  
  times = runif(n_sub_cnas, t_subclone_mrca, t_final) %>% sort()
  
  subclonal_cnas = data.frame()
  for (i in 1:n_sub_cnas){
    valid_k = FALSE
    while (valid_k == FALSE){
      karyotypes = sample_karyo()
      #print(karyotypes)
      model = sample_model(karyotypes, purity, ccf) #%>% pull(model_id) %>% unique()
      if (typeof(model) != 'character'){
        model = model %>% pull(model_id) %>% unique()
        valid_k = TRUE
        }
      }
    
    chrs = sample_chr(filtered_chr = chromosomes_without_cnas)
    chromosomes_without_cnas = chromosomes_without_cnas[!(chromosomes_without_cnas==chrs)]
    form_to = sample_from_to(chrs)
    from = form_to[1]
    to = form_to[2]
    
    sc = data.frame('chr' = chrs, 'from'=from, 'to'= to, 'karyotype'= karyotypes, 'model_id'= model)
    #print(sc)
    subclonal_cnas = rbind(sc, subclonal_cnas)
  }
  subclonal_cnas = subclonal_cnas %>% mutate('times'= times) %>% rowwise() %>% 
    mutate('Major'=1,'minor'=1,'Major2'= strsplit(karyotype, ':')[[1]][1], 'minor2' = strsplit(karyotype, ':')[[1]][2]) %>%
    mutate(segment_id = paste0(chr,':' ,from,':' ,to, ':',Major, ':',minor, ':', ccf), CCF = ccf)
  subclonal_cnas
}

merge_clonal_subclonal_cnas = function(cnas, sub_cnas){
  
  subclonal_chrs = sub_cnas$chr
  sub_cnas = sub_cnas %>% mutate(genotype= model_id)
  
  new_rows = data.frame()
  for (i in 1:length(subclonal_chrs)){
    new_row_11 = data.frame('chr' = subclonal_chrs[i], 'from'=0, 
                            'to'= sub_cnas[i,]$from, 'karyotype'= '1:1', 'genotype'= 'A1B1', 
                            'times'= 0, model_id='clonal', Major2='na',minor2='na',  Major=1,minor=1, 'CCF'=1) %>%
      mutate(segment_id = paste0(chr,':' ,from,':' ,to, ':',Major, ':',minor, ':1'))
    new_row_12 = data.frame('chr' = subclonal_chrs[i], 'from'=sub_cnas[i,]$to, 'to'= get_chr_length(subclonal_chrs[i]), 
                            'karyotype'= '1:1', 'genotype'= 'A1B1', 'times'= 0, 
                            model_id='clonal', Major2='na',minor2='na',Major=1,minor=1,'CCF'=1) %>%
      mutate(segment_id = paste0(chr,':' ,from,':' ,to, ':',Major, ':',minor, ':1'))
    new_rows = rbind(new_rows, new_row_11, new_row_12)
  }
  
  cnas = cnas %>% mutate('Major2'= 'na', 'minor2'='na', 'model_id'= 'clonal', 'CCF'=1) %>% filter(!(chr %in% subclonal_chrs))
  merged_cnas = rbind(new_rows, cnas,sub_cnas)
  merged_cnas
}

# 4. SIMULATE SUBCLONAL SNVs :
simulate_subclonal_snvs_df = function(sub_cnas, mu= 1e-8, w= 2, t_f = 10, coverage = 100, purity=1){
  snvs_df = data.frame()
  for (i in 1:length(sub_cnas$chr)){
    n_muts = rpois(1, 2*mu*w*t_f*get_chr_length(sub_cnas$chr[i]))
    kc = paste0('1:1-', sub_cnas$karyotype[i])
    peaks = get_models(karyotype_combination= kc, purity= purity, ccf= sub_cnas$CCF[i]) %>% 
      filter(model_id==sub_cnas$model_id[i]) %>% pull(peak)
    dp = rpois(n_muts, coverage)
    nv = rbinom(n_muts, dp, peaks)
    from = sample.int(sub_cnas$to[i] - sub_cnas$from[i], n_muts)
    to = from+1
    
    snvs_df_i = data.frame('chr' = sub_cnas$chr[i], 'from'= from, 'to'= to, 
                           'Major'= sub_cnas$Major2[i],'minor'=sub_cnas$minor2[i],'NV'=nv,'DP'=dp, 'VAF'=nv/dp,
                           'allele'= '', 'segment_id'= sub_cnas$segment_id[i])
    snvs_df= rbind(snvs_df_i, snvs_df)
  }
  snvs_df
}

# 5. ADD TAIL MUTATIONS
simulate_tail_snvs = function(all_cnas, mu= 1e-8, w= 2, t_f = 10, coverage = 100, purity=1){

  chromosomes = c(paste0('chr', 1:22), 'chrX', 'chrY')
  len = sapply(chromosomes, function(x){
    get_chr_length(x)
  })
  len = sum(len)

  segment_ids = all_cnas$segment_id #%>% pull(segment_id)

  tail_prop =.5

  n_muts_tail_1 <- rpois(1, len * mu * w * t_f * tail_prop)
  n_muts_tail_2 <- rpois(1, len * mu * w * t_f * tail_prop)

  dp_tail_1 <- rpois(n_muts_tail_1, coverage)
  dp_tail_2 <- rpois(n_muts_tail_2, coverage)

  tail_nv_1 <- round(my_pareto_betabin(n_muts_tail_1) * dp_tail_1, 1)
  tail_nv_2 <- round(my_pareto_betabin(n_muts_tail_2) * dp_tail_2, 1)
  
  require(tidyverse)

  tail_1 <- data.frame('chr' = '', 'from'= '', 'to'= '',
                       'Major'= 1,'minor'= 1,'NV'=tail_nv_1,'DP'=dp_tail_1, 'VAF'=tail_nv_1/dp_tail_1,
                       'allele'= '', 'segment_id'= segment_ids[sample.int(length(segment_ids), n_muts_tail_1, replace = TRUE)]) %>%
    rowwise %>% mutate('chr'= strsplit(segment_id, ':')[[1]][1],
                        'Major'= as.integer(strsplit(segment_id, ':')[[1]][4]),
                        'minor'= as.integer(strsplit(segment_id, ':')[[1]][5]))

  tail_2 <- data.frame('chr' = '', 'from'= '', 'to'= '',
                       'Major'= 1,'minor'= 1,'NV'=tail_nv_2,'DP'=dp_tail_2, 'VAF'=tail_nv_2/dp_tail_2,
                       'allele'= '', 'segment_id'= segment_ids[sample.int(length(segment_ids), n_muts_tail_2, replace = TRUE)]) %>%
    rowwise %>% mutate('chr'= strsplit(segment_id, ':')[[1]][1],
                        'Major'= as.integer(strsplit(segment_id, ':')[[1]][4]),
                        'minor'= as.integer(strsplit(segment_id, ':')[[1]][5]))


  return(rbind(tail_1, tail_2))


}


# 6. ADD SNVS minding the allele composition
simulate_snps_seg = function(segment, bin_size= 50000, purity= 1, ploidy=2){
  n_bins <- (segment$to - segment$from) / bin_size
  n_baf = rpois(n_bins, (1/1000)*bin_size)
  
  if (segment$CCF==1){
    exp_baf <- clonal_expected_baf(segment$karyotype, purity)
    exp_dr <- clonal_expected_dr(segment$karyotype, purity,ploidy)
  }else{
    kc= paste0('1:1-', segment$karyotype)
    g1 = get_models(karyotype_combination= kc, purity= purity, ccf=segment$CCF) %>% 
      filter(model_id==segment$model_id) %>% pull(genotype_1)
    g2 = get_models(karyotype_combination= kc, purity= purity, ccf=segment$CCF) %>% 
      filter(model_id==segment$model_id) %>% pull(genotype_2)
    exp_baf= expected_baf(k1='1:1', k2=segment$karyotype, purity=purity, ccf=segment$CCF, g1, g2)
    exp_dr= expected_dr(k1='1:1', k2=segment$karyotype, purity=purity, ccf=segment$CCF, ploidy = ploidy)
  }
  

  alpha <- ((n_baf - 2) * exp_baf + 1) / (1 - exp_baf)
  baf <- rbeta(n_bins, shape1 = alpha, shape2 = n_baf) 
  dr <- rgamma(n_bins, shape = exp_dr * sqrt(n_baf) + 1, rate = sqrt(n_baf))
  
  SNP_df <- data.frame(
    chr = segment$chr,
    from = seq(from=1, to=get_chr_length(segment$chr), by=50000)[1:n_bins],
    to = seq(from=50000, to=get_chr_length(segment$chr), by=50000)[1:n_bins],
    BAF = baf, N.BAF = n_baf, DR = dr
  )

  return(SNP_df)
}

simulate_snps = function(all_cnas, bin_size= 50000, purity= 1, ploidy=2){
  all_snps = lapply(all_cnas$segment_id, function(x){
    current_cna = all_cnas %>% filter(segment_id==x)
    snps = simulate_snps_seg(current_cna, bin_size= bin_size, purity= purity, ploidy=ploidy)
    snps
  })
  all_snps = Reduce(rbind, all_snps)
  all_snps
}

# 7. GET THE CNAQC SAMPLE

simulate_sample = function(t_subclone_mrca = 5, cna_rate = 5,
                           mu= 1e-8, w= 2, coverage = 100, purity=1,
                           sub_cnas_rate= 3, t_final= 10, ccf = .5, bin_size= 50000){
  cnas = simulate_clonal_cnas(t_subclone_mrca = t_subclone_mrca, cna_rate = cna_rate)
  cl_snvs = simulate_clonal_snvs_df(cnas=cnas,
                                    mu= mu, w= w, t_subclone_mrca = t_subclone_mrca, 
                                    coverage = coverage, purity=purity)
  sub_cnas = simulate_subclonal_cnas(cna=cnas,
                                     sub_cnas_rate= sub_cnas_rate, t_subclone_mrca = t_subclone_mrca, 
                                     t_final= t_final, ccf = ccf, purity =purity)
  all_cnas = merge_clonal_subclonal_cnas(cnas, sub_cnas)
  subclonal_muts = simulate_subclonal_snvs_df(sub_cnas, mu= mu, w= w, t_f = t_f, coverage = coverage, purity=purity)
  tail_muts = simulate_tail_snvs(all_cnas, mu= mu, w= w, t_f = t_f, coverage = coverage, purity= purity)
  all_snvs = rbind(cl_snvs, subclonal_muts, tail_muts)
  ploidy= compute_approx_ploidy(all_cnas)
  all_snps = simulate_snps(all_cnas, bin_size= bin_size, purity= purity, ploidy=ploidy)
  
  return(list('cna'= all_cnas, 'mutations'= all_snvs, 'snps'= all_snps, 'purity'= purity))
}


#hist(all_snvs %>% filter(Major==2, minor==1) %>% pull(VAF), breaks =50)



