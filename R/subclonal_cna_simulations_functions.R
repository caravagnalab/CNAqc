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
sample_bases = function(){
  b = c('A', 'C', 'T', 'G')
  b[sample.int(4, 2)]
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
  from = as.integer(sample.int(chr_l-l, 1))
  to = as.integer(from + l)
  
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

# simulate_clonal_cnas = function(t_subclone_mrca = 5, cna_rate = 5){
#   n_cnas = rpois(1, cna_rate)
#   times = runif(n_cnas, 0, t_subclone_mrca) %>% sort()
#   
#   clonal_cnas = data.frame()
#   for (i in 1:n_cnas){
#     karyotypes = sample_karyo()
#     chrs = sample_chr()
#     form_to = sample_from_to(chrs)
#     from = form_to[1]
#     to = form_to[2]
#     genotype = sample_genotype(karyotypes) 
#     cc = data.frame('chr' = chrs, 'from'=from, 'to'= to, 'karyotype'= karyotypes, 'genotype'= genotype)
#     clonal_cnas = rbind(cc, clonal_cnas)
#   }
#   clonal_cnas = clonal_cnas %>% mutate('times'= times) %>% fill_the_gaps() %>% rowwise() %>% 
#     mutate('Major'= strsplit(karyotype, ':')[[1]][1], 'minor' = strsplit(karyotype, ':')[[1]][2]) %>%
#     mutate(segment_id = paste0(chr,':' ,from,':' ,to, ':',Major, ':',minor, ':1'))
# }

simulate_clonal_cnas = function(t_subclone_mrca = 5, cna_rate = 5){
  n_cnas = rpois(1, cna_rate)
  times = runif(n_cnas, 0, t_subclone_mrca) %>% sort()

  clonal_cnas = get_reference(ref= 'hg19') %>% select(chr, length) %>% 
    mutate(from=0, to= length, karyotype= '1:1', genotype = 'A1B1', times =0, Major=1, minor=1) %>% select(!length)
  
  for (i in 1:n_cnas){
    karyotypes = sample_karyo()
    chrs = sample_chr()
    form_to = sample_from_to(chrs)
    from_i = form_to[1]
    to_i = form_to[2]
    genotype = sample_genotype(karyotypes) 
    clonal_cnas = clonal_cnas %>% filter(chr!=chrs)
    
    if (from_i ==0){
      r1 = data.frame(chr= chrs, from=from_i, to= to_i, karyotype= karyotypes, genotype = genotype, 
                      times = times[i], Major=strsplit(karyotypes, ':')[[1]][1], minor=strsplit(karyotypes, ':')[[1]][2])
      r2 = data.frame(chr= chrs, from=to_i+1, to= get_chr_length(chrs), karyotype= '1:1', genotype = 'A1B1', times =0, Major=1, minor=1)
      r = rbind(r1,r2)
    }
    if (to_i == get_chr_length(chrs)){
      r1 = data.frame(chr= chrs, from=0, to= from_i-1, karyotype= '1:1', genotype = 'A1B1', times =0, Major=1, minor=1)
      r2 = data.frame(chr= chrs, from=from_i, to= to_i, karyotype= karyotypes, genotype = genotype, 
                      times = times[i], Major=strsplit(karyotypes, ':')[[1]][1], minor=strsplit(karyotypes, ':')[[1]][2])
      r = rbind(r1,r2)
    }
    if (from_i > 0 & to_i < get_chr_length(chrs)){
      r1 = data.frame(chr= chrs, from=0, to= from_i-1, karyotype= '1:1', genotype = 'A1B1', times =0, Major=1, minor=1)
      r2 = data.frame(chr= chrs, from=from_i, to= to_i, karyotype= karyotypes, genotype = genotype, 
                      times = times[i], Major=strsplit(karyotypes, ':')[[1]][1], minor=strsplit(karyotypes, ':')[[1]][2])
      r3 = data.frame(chr= chrs, from= to_i+1, to= get_chr_length(chrs), karyotype= '1:1', genotype = 'A1B1', times =0, Major=1, minor=1)
      r = rbind(r1, r2,r3)
    }
    
    
    clonal_cnas = rbind(clonal_cnas, r)
  }
  
  clonal_cnas= clonal_cnas #%>% mutate(segment_id = paste0(chr,':' ,from,':' ,to, ':',Major, ':',minor, ':1'))
}

# 2. SIMULATE CLONAL SNVs 
# 2.1 Place the clonal snvs on each segment :
#         - the older ones with higher multiplicities will accumulate with rate depending on 
#           mutation rate and growth rate, over a period of time going from 0 to time of segment formation
#         - the newer ones with rate depending on a period going from the time of segment formation to subclone birth

generate_clonal_snvs = function(cna, seg_id, karyotype, chromosome, #al='A1B1',
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
  #p = rbeta(n_snvs, round(100*peak_mean), round(100*(1-peak_mean)))
  p = peak_mean
  nv = rbinom(n_snvs, dp, p)
  from_i = runif(n_snvs, cna$from, cna$to-1)
  #from_i = sample.int(get_chr_length(chromosome)-1, n_snvs)
  to_i = from_i +1
  #allele = c()
  ref= c()
  alt = c()
  for (i in 1:n_snvs){
    bases = sample_bases()
    ref[i] = bases[1]
    alt[i] = bases[2]
  }
  #for (i in 1:n_snvs){allele[i]=sample_allele(al)}
  #allele = rep(sample_allele('A1B1'), n_snvs)
  maj = strsplit(karyotype, ':')[[1]][1]
  minor = strsplit(karyotype, ':')[[1]][2]
  snvs_df = data.frame('chr' = chromosome, 'from'=from_i, 'to'= to_i, 'ref'= ref, 'alt'=alt,
                       'Major'= maj,'minor'=minor,'NV'=nv,'DP'=dp, 'VAF'=nv/dp,
                       #'allele'= allele, 
                       'segment_id'= seg_id)
  snvs_df
  
}

simulate_clonal_snvs_df = function(cnas, mu= 1e-8, w= 2, t_subclone_mrca = 5, coverage = 100, purity=1){
  segments = cnas %>% filter(model_id == 'clonal') %>% pull(segment_id)
  clonal_snvs_df = lapply(segments, function(x){
    chromosome = strsplit(x, ':')[[1]][1]
    snvs_df = data.frame()
    cna_s = cnas %>% filter(segment_id ==x)
    k = paste0(cna_s$Major, ':', cna_s$minor)
    #muts=0
    #while (muts < 50){
      #print(x)
    snvs_df_i = generate_clonal_snvs(cna = cnas %>% filter(segment_id == x), seg_id=x, karyotype= k,
                                     chromosome=chromosome, #al= cna_s$genotype,
                                     mu=mu, w= w, t_subclone_mrca = t_subclone_mrca,
                                     coverage = coverage, purity=purity, t_cna=cna_s$times)
    
    #muts = length(snvs_df_i$chr)
    #}
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
    mutate('Major'=1,'minor'=1,'Major_2'= strsplit(karyotype, ':')[[1]][1], 'minor_2' = strsplit(karyotype, ':')[[1]][2], CCF = ccf) #%>%
    #mutate(segment_id = paste0(chr,':' ,from,':' ,to, ':',Major, ':',minor, ':', ccf), CCF = ccf)
  subclonal_cnas
}

merge_clonal_subclonal_cnas = function(cnas, sub_cnas){
  
  subclonal_chr = sub_cnas %>% pull(chr)
  cnas = cnas %>% filter(!(chr %in% subclonal_chr))
  
  new_subclonal_cnas = data.frame()
  for (i in 1:length(subclonal_chr)){
    from_i= sub_cnas[i,]$from
    to_i= sub_cnas[i,]$to
    chrs= sub_cnas[i,]$chr
    if (from_i ==0){
      r1 = sub_cnas[i,]
      r2 = sub_cnas[i,] %>% mutate(from = to_i+1, to = get_chr_length(sub_cnas[i,]$chr), karyotype='1:1', 
                                   Major_2=0, minor_2=0, CCF=1, times=0)
      r = rbind(r1,r2)
    }
    if (to_i == get_chr_length(chrs)){
      r1 = sub_cnas[i,] %>% mutate(to= from_i-1) %>% mutate(from = 0, karyotype='1:1', 
                                   Major_2=0, minor_2=0, CCF=1, times=0)
      r2 = sub_cnas[i,]
      r = rbind(r1,r2)
    }
    if (from_i > 0 & to_i < get_chr_length(chrs)){
      r1 = sub_cnas[i,] %>% mutate(to =from_i-1) %>% mutate(from = 0, karyotype='1:1', 
                                   Major_2=0, minor_2=0, CCF=1, times=0)
      r2 = sub_cnas[i,]
      r3 = sub_cnas[i,] %>% mutate(from = to_i+1, to = get_chr_length(sub_cnas[i,]$chr), karyotype='1:1', 
                                   Major_2=0, minor_2=0, CCF=1, times=0)
      r = rbind(r1, r2,r3)
    }
    new_subclonal_cnas = rbind(new_subclonal_cnas, r)
  }
  cnas= cnas %>% mutate(model_id='clonal', Major_2=0, minor_2=0, CCF=1) %>% select(!genotype)
  all_cnas = rbind(cnas, new_subclonal_cnas)
}

# 4. SIMULATE SUBCLONAL SNVs :
simulate_subclonal_snvs_df = function(sub_cnas, mu= 1e-8, w= 2, t_f = 10, coverage = 100, purity=1){
  sub_cnas = sub_cnas %>% filter(CCF<1)
  snvs_df = data.frame()
  #sub_cnas_ids
  for (i in 1:length(sub_cnas$chr)){
    n_muts = rpois(1, 2*mu*w*t_f*get_chr_length(sub_cnas$chr[i]))
    k = paste0(sub_cnas$Major_2[i], ':', sub_cnas$minor_2[i])
    kc = paste0('1:1-', k)
    peaks = get_models(karyotype_combination= kc, purity= purity, ccf= sub_cnas$CCF[i]) %>% 
      filter(model_id==sub_cnas$model_id[i]) %>% pull(peak)
    dp = rpois(n_muts, coverage)
    nv = rbinom(n_muts, dp, peaks)
    from_i = runif(n_muts, sub_cnas$from[i], sub_cnas$to[i]-1 ) #sample.int(sub_cnas$to[i] - sub_cnas$from[i], n_muts)
    to_i = from_i+1
    
    ref= c()
    alt = c()
    for (j in 1:n_muts){
      bases = sample_bases()
      ref[j] = bases[1]
      alt[j] = bases[2]
    }
    
    c = sub_cnas$chr[i]
    snvs_df_i = data.frame('chr' = sub_cnas$chr[i], 'from'= from_i, 'to'= to_i, 'ref' = ref, 'alt'= alt,
                           #'Major'= sub_cnas$Major_2[i],'minor'=sub_cnas$minor_2[i],
                           'Major'= 1,'minor'=1,
                           'NV'=nv,'DP'=dp, 'VAF'=nv/dp,
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
  
  ref= c()
  alt = c()
  for (i in 1:n_muts_tail_1){
    bases = sample_bases()
    ref[i] = bases[1]
    alt[i] = bases[2]
  }
  
  chrms = c()
  from = c()
  to = c()
  for (j in 1:n_muts_tail_1){
    chrms[j] = sample_chr()
    from[j]= runif(1, 0, get_chr_length(chrms[j]))
    to[j] = from[j] +1
  }
  
  require(tidyverse)

  tail_1 <- data.frame('chr' = chrms, 'from'= from, 'to'= to, 'ref'= ref, 'alt'=alt, 
                       'Major'= 1,'minor'= 1,'NV'=tail_nv_1,'DP'=dp_tail_1, 'VAF'=tail_nv_1/dp_tail_1,
                       'allele'= '', 'segment_id'= segment_ids[sample.int(length(segment_ids), n_muts_tail_1, replace = TRUE)]) %>%
    rowwise %>% mutate('chr'= strsplit(segment_id, ':')[[1]][1],
                        'Major'= as.integer(strsplit(segment_id, ':')[[1]][4]),
                        'minor'= as.integer(strsplit(segment_id, ':')[[1]][5]))
  
  ref= c()
  alt = c()
  for (i in 1:n_muts_tail_2){
    bases = sample_bases()
    ref[i] = bases[1]
    alt[i] = bases[2]
  }
  
  chrms = c()
  from = c()
  to = c()
  for (j in 1:n_muts_tail_2){
    chrms[j] = sample_chr()
    from[j]= runif(1, 0, get_chr_length(chrms[j]))
    to[j] = from[j] +1
  }

  tail_2 <- data.frame('chr' = chrms, 'from'= from, 'to'= to, 'ref'= ref, 'alt'=alt, 
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
  
  if (segment$model_id=='clonal'){
    exp_baf <- clonal_expected_baf(paste0(segment$Major, ':', segment$minor), purity)
    exp_dr <- clonal_expected_dr(paste0(segment$Major, ':', segment$minor), purity,ploidy)
  }else{
    kc= paste0('1:1-', segment$Major_2, ':', segment$minor_2)
    g1 = get_models(karyotype_combination= kc, purity= purity, ccf=segment$CCF) %>% 
      filter(model_id==segment$model_id) %>% filter(!is.na(genotype_1)) %>% pull(genotype_1) %>% unique()
    g2 = get_models(karyotype_combination= kc, purity= purity, ccf=segment$CCF) %>% 
      filter(model_id==segment$model_id) %>% filter(!is.na(genotype_2)) %>% pull(genotype_2) %>% unique()
    exp_baf= min(expected_baf(k1='1:1', paste0(segment$Major_2, ':', segment$minor_2), purity=purity, ccf=segment$CCF, g1, g2))
    
    exp_dr= expected_dr(k1='1:1', k2= paste0(segment$Major_2, ':', segment$minor_2), purity=purity, ccf=segment$CCF, ploidy = ploidy)
  }
  

  alpha <- ((n_baf - 2) * exp_baf + 1) / (1 - exp_baf)
  baf <- rbeta(n_bins, shape1 = alpha, shape2 = n_baf) 
  dr <- rgamma(n_bins, shape = exp_dr * sqrt(n_baf) + 1, rate = sqrt(n_baf))
  
  SNP_df <- data.frame(
    chr = segment$chr,
    from = seq(from=segment$from, to=segment$to - 50000, by=50000)[1:n_bins],
    to = seq(from=segment$from+50000, to=segment$to, by=50000)[1:n_bins],
    BAF = baf, N.BAF = n_baf, DR = dr,
    segment_id = segment$segment_id
  )

  return(SNP_df)
}


simulate_snps = function(all_cnas, bin_size= 50000, purity= 1, ploidy=2){
  #chromosomes = c(paste0('chr', 1:22), 'chrX', 'chrY')
  all_snps = lapply(all_cnas$segment_id, function(x){
    current_cna = all_cnas %>% filter(segment_id==x)
    snps = simulate_snps_seg(current_cna, bin_size= bin_size, purity= purity, ploidy=ploidy)
    snps
  })
  
  all_snps = Reduce(rbind, all_snps)
  all_snps
}

# 7. GET THE CNAQC SAMPLE

# t_subclone_mrca = 5
# cna_rate = 5
# mu= 1e-8
# w= 2
# coverage = 100
# purity=1
# sub_cnas_rate= 3
# t_f= 10
# ccf = .5
# bin_size= 50000

simulate_sample = function(t_subclone_mrca = 5, cna_rate = 5,
                           mu= 1e-8, w= 2, coverage = 100, purity=1,
                           sub_cnas_rate= 3, t_f= 10, ccf = .5, bin_size= 50000){
  cnas = simulate_clonal_cnas(t_subclone_mrca = t_subclone_mrca, cna_rate = cna_rate)
  sub_cnas = simulate_subclonal_cnas(cna=cnas,
                                     sub_cnas_rate= sub_cnas_rate, t_subclone_mrca = t_subclone_mrca, 
                                     t_final= t_f, ccf = ccf, purity =purity)
  all_cnas = merge_clonal_subclonal_cnas(cnas, sub_cnas) %>% select(!karyotype) %>% mutate(segment_id= paste0(chr,':' ,from,':' ,to, ':',Major, ':',minor, ':', CCF))
  
  cl_snvs = simulate_clonal_snvs_df(cnas=all_cnas,
                                    mu= mu, w= w, t_subclone_mrca = t_subclone_mrca, 
                                    coverage = coverage, purity=purity) 
  subclonal_muts = simulate_subclonal_snvs_df(all_cnas, mu= mu, w= w, t_f = t_f, coverage = coverage, purity=purity) %>% select(!allele)
  tail_muts = simulate_tail_snvs(all_cnas, mu= mu, w= w, t_f = t_f, coverage = coverage, purity= purity) %>% select(!allele)
  all_snvs = rbind(cl_snvs, subclonal_muts, tail_muts) %>% mutate(NV = as.integer(NV), DP = as.integer(DP)) 
  ploidy= compute_approx_ploidy(all_cnas)
  all_cnas = all_cnas %>% mutate(model_id = ifelse(CCF==1, 'clonal', model_id))
  all_snps = simulate_snps(all_cnas, bin_size= bin_size, purity= purity, ploidy=ploidy)
  
  return(list('cna'= all_cnas, 'mutations'= all_snvs, 'snps'= all_snps, 'purity'= purity, 'reference_genome'='hg19'))
}


# 8. CHECK SOLUTION

# Expected BAF clonal and subclonal distributions
clonal_baf_ditribution = function(n, k1, purity, n_points) {
  E_baf <- clonal_expected_baf(k1, purity)
  alpha <- ((n - 2) * E_baf + 1) / (1 - E_baf)
  s = rbeta(n_points, shape1 = alpha, shape2 = n)
  #s <- dbeta(baf_obs, shape1 = alpha, shape2 = n) 
  return(s)
}
subclonal_baf_distribution = function(n, k1, k2, purity, ccf, g1, g2, n_points){
  E_baf <- expected_baf(k1, k2, purity, ccf, g1, g2)
  alpha <- ((n - 2) * E_baf + 1) / (1 - E_baf)
  s <- rbeta(n_points, shape1 = alpha, shape2 = n) 
  #s <- dbeta(baf_obs, shape1 = alpha, shape2 = n) 
  return(s)
}
expected_baf_ditribution= function(n, k1, k2='', purity, ccf=1, g1='', g2='', n_points){
  if (ccf ==1){
    expected_baf_distr = clonal_baf_ditribution(n, k1, purity, n_points)
  }else{
    expected_baf_distr = subclonal_baf_distribution(n, k1, k2, purity, ccf, g1, g2, n_points)
  }
  expected_baf_distr
}
# Plot Baf expecteed vs observed distributions
plot_baf_exp_vs_obs = function(patched_x, all_res_seg, n_snps, n_points, which, s= 5){

  # compute expected baf distribution according to model (given, correct, best)
  if (which=='given'){
    sim_color='steelblue'
    expected_baf_d = expected_baf_ditribution(n= n_snps, k1= all_res_seg$k1[1], k2='', 
                                            purity =  patched_x$purity, ccf=1, g1='', g2='', n_points)
  }
  
  if (which == 'correct'){
    sim_color='forestgreen'
    if (all_res_seg$k2[1]==''){
      expected_baf_d = expected_baf_ditribution(n= n_snps, k1= all_res_seg$k1[1], k2='', 
                                                purity =  patched_x$purity, ccf=1, g1='', g2='', n_points)
    }else{
      genotypes = get_models(karyotype_combination = paste0(all_res_seg$k1[1],'-' ,all_res_seg$k2[1]),
                             purity = x$purity, ccf = all_res_seg$CCF_sim[1]) %>% filter(model_id == all_res_seg$model_id[1]) %>% 
        select(genotype_1,genotype_2) 
      g1 = genotypes %>% filter(!is.na(genotype_1)) %>% pull(genotype_1) %>% unique()
      g2 = genotypes %>% filter(!is.na(genotype_2)) %>% pull(genotype_2) %>% unique()
      expected_baf_d = expected_baf_ditribution(n= n_snps, k1= all_res_seg$k1[1], k2= all_res_seg$k2[1], 
                                                purity =  patched_x$purity, ccf= all_res_seg$CCF_sim[1], g1, g2, n_points)
    }
    
  }
  
  if (which == 'inferred'){
    sim_color='indianred'
    if (all_res_seg$k2_best[1]==''){
      expected_baf_d = expected_baf_ditribution(n= n_snps, k1= all_res_seg$k1_best[1], k2='', 
                                                purity =  patched_x$purity, ccf=1, g1='', g2='', n_points)
    }else{
      genotypes = get_models(karyotype_combination = paste0(all_res_seg$k1_best[1],'-' ,all_res_seg$k2_best[1]),
                             purity = x$purity, ccf = all_res_seg$CCF_best[1]) %>% filter(model_id == all_res_seg$model_id_best[1]) %>% 
        select(genotype_1,genotype_2) 
      g1 = genotypes %>% filter(!is.na(genotype_1)) %>% pull(genotype_1) %>% unique()
      g2 = genotypes %>% filter(!is.na(genotype_2)) %>% pull(genotype_2) %>% unique()
      expected_baf_d = expected_baf_ditribution(n= n_snps, k1= all_res_seg$k1_best[1], k2= all_res_seg$k2_best[1], 
                                                purity =  patched_x$purity, ccf= all_res_seg$CCF_best[1], g1, g2, n_points)
    }
    
  }
  
  e_baf_d = patched_x$snps %>% filter(segment_id==all_res_seg$segment_id[1]) %>% 
    select(from, to) %>% mutate(exp_baf = expected_baf_d) %>% as_tibble() %>% mutate(from=as.integer(from), to=as.integer(to))
  
  # the observed baf distribution is constant 
  observed_baf_d = patched_x$snps %>% filter(segment_id==all_res_seg$segment_id[1]) %>% select(from, to, BAF)
  
  ggplot() +
    geom_segment(data= e_baf_d, aes(x = from, xend= to, y= exp_baf, yend = exp_baf), color= sim_color, size = s)+
    #geom_density_2d(data= e_baf_d, aes(x = from, y= exp_baf))+
    geom_segment(data = observed_baf_d %>% as_tibble(), aes(x = from, xend= to, y= BAF, yend = BAF), color= 'black', size = s)+
    theme_bw() + labs(x= '', y= 'BAF') + ylim(0,1)
}

# Expected DR clonal and subclonal distributions
clonal_dr_distribution = function(n, k1, purity, ploidy=2, n_points){
  expected_dpr <- clonal_expected_dr(k1, purity, ploidy)
  dr <- rgamma(n_points, shape = expected_dpr * sqrt(n) + 1, rate = sqrt(n))
  return(dr)
}
subclonal_dr_distribution = function(n, k1, k2, purity, ccf, ploidy = 2, n_points){
  expected_dpr <- expected_dr(k1, k2, purity, ccf, ploidy)
  dr <- rgamma(n_points, shape = expected_dpr * sqrt(n) + 1, rate = sqrt(n))
  return(dr)
}
expected_dr_distribution = function(n, k1, k2='', purity, ccf=1, ploidy = 2, n_points){
  if (ccf==1){
    distr = clonal_dr_distribution(n, k1, purity, ploidy, n_points)
  }else{
    distr = subclonal_dr_distribution(n, k1, k2, purity, ccf, ploidy, n_points)
  }
}
# Plot dr expecteed vs observed distributions
plot_dr_exp_vs_obs = function(patched_x, all_res_seg, n_snps, n_points, which, s= 5){
  
  # compute expected baf distribution according to model (given, correct, best)
  if (which=='given'){
    sim_color = 'steelblue'
    expected_dr_d = expected_dr_distribution(n= n_snps, k1= all_res_seg$k1[1], k2='', 
                                              purity =  patched_x$purity, ccf=1, ploidy = patched_x$ploidy, n_points)
  }
  
  if (which == 'correct'){
    sim_color = 'forestgreen'
    if (all_res_seg$k2[1]==''){
    expected_dr_d = expected_dr_distribution(n= n_snps, k1= all_res_seg$k1[1], k2='', 
                                             purity =  patched_x$purity, ccf=1, ploidy = patched_x$ploidy, n_points)
    }else{
      expected_dr_d = expected_dr_distribution(n= n_snps, k1= all_res_seg$k1[1], k2=all_res_seg$k2[1], 
                                               purity =  patched_x$purity, ccf=all_res_seg$CCF_sim[1], 
                                               ploidy = patched_x$ploidy, n_points)
    }
  }
  
  if (which == 'inferred'){
    sim_color = 'indianred'
    if (all_res_seg$k2_best[1]==''){
      expected_dr_d = expected_dr_distribution(n= n_snps, k1= all_res_seg$k1_best[1], k2='', 
                                               purity =  patched_x$purity, ccf=1, ploidy = patched_x$ploidy, n_points)
    }else{
      expected_dr_d = expected_dr_distribution(n= n_snps, k1= all_res_seg$k1_best[1], k2=all_res_seg$k2_best[1], 
                                               purity =  patched_x$purity, ccf=all_res_seg$CCF_best[1], 
                                               ploidy = patched_x$ploidy, n_points)
    }
    
  }
  
  e_dr_d = patched_x$snps %>% filter(segment_id==all_res_seg$segment_id[1]) %>% 
    select(from, to) %>% mutate(exp_dr = expected_dr_d) %>% as_tibble() %>% mutate(from=as.integer(from), to=as.integer(to))
  # the observed baf distribution is constant 
  observed_dr_d = patched_x$snps %>% filter(segment_id==all_res_seg$segment_id[1]) %>% select(from, to, DR)
  
  ggplot() +
    geom_segment(data= e_dr_d, aes(x = from, xend= to, y= exp_dr, yend = exp_dr), color= sim_color, size=s)+
    geom_segment(data = observed_dr_d %>% as_tibble(), aes(x = from, xend= to, y= DR, yend = DR), color= 'black',size=s)+
    theme_bw() + labs(x= '', y= 'DR')
}

# Expected VAF distribution
vaf_distribution = function(n_points, coverage = NULL, peaks){
  dp = c()
  nv = c()
  for (i in 1:length(peaks)){
    dp = c(dp, rpois(as.integer(n_points/length(peaks)), coverage))
    nv = c(nv, rbinom(as.integer(n_points/length(peaks)), size=dp, prob=peaks[i]))
  }
  
  return(data.frame('nv'= nv, 'dp'= dp))
}
plot_vaf_exp_vs_obs = function(patched_x, all_res_seg, n_points, coverage, which){
  if (which == 'given'){
    sim_color ='steelblue'
    peaks = get_clonal_peaks(all_res_seg$k1[1], patched_x$purity)
  }
  if (which == 'correct'){
    sim_color='forestgreen'
    if (all_res_seg$k2[1]==''){
      peaks = get_clonal_peaks(all_res_seg$k1[1], patched_x$purity)
    }else{
    peaks = get_models(karyotype_combination = paste0(all_res_seg$k1[1],'-' ,all_res_seg$k2[1]),
                       purity = x$purity, ccf = all_res_seg$CCF_sim[1]) %>% filter(model_id == all_res_seg$model_id[1]) %>% pull(peak)
    }
  }
  if (which == 'inferred'){
    sim_color='indianred'
    peaks = all_res_seg %>% pull(peak)
  }
  vd = vaf_distribution(n_points, coverage, peaks)
  ggplot(vd) +
    geom_histogram(aes(x= nv/dp), fill= sim_color, alpha=.5)+
    geom_histogram(data = patched_x$mutations %>% filter(segment_id==all_res_seg$segment_id), aes(x= VAF), fill= 'black', alpha=.5)+
    theme_bw()+labs(y='', x='VAF') + xlim(0,1)
}

# Get likelihood Given and Correct models 
get_likelihood_simulated_model = function(patched_x, all_res_seg, muts, snps, which){
  
  if (which == 'given'){
  # Given solution LL
  v <- VAF_LL(NV= muts$NV, DP= muts$DP, peaks = data.frame('peak'=get_clonal_peaks(all_res_seg$k1[1], purity=patched_x$purity)))
  b = sum(log(clonal_baf_ll(baf_obs= snps$BAF, n = snps$N.BAF, k1=all_res_seg$k1[1], purity=patched_x$purity)))
  d <- sum(log(clonal_dr_ll(snps$DR, n = snps$N.BAF, all_res_seg$k1[1], patched_x$purity, patched_x$ploidy)))
}else{
  if (all_res_seg$k2[1]==''){
    v <- VAF_LL(NV= muts$NV, DP= muts$DP, peaks = data.frame('peak'=get_clonal_peaks(all_res_seg$k1[1], purity=patched_x$purity)))
    b = sum(log(clonal_baf_ll(baf_obs= snps$BAF, n = snps$N.BAF, k1=all_res_seg$k1[1], purity=patched_x$purity)))
    d <- sum(log(clonal_dr_ll(snps$DR, n = snps$N.BAF, all_res_seg$k1[1], patched_x$purity, patched_x$ploidy)))
  }else{
    peaks = data.frame('peak' = get_models(karyotype_combination = paste0(all_res_seg$k1[1],'-' ,all_res_seg$k2[1]),
                       purity = x$purity, ccf = all_res_seg$CCF_sim[1]) %>% filter(model_id == all_res_seg$model_id[1]) %>% pull(peak))
    v <- VAF_LL(NV = muts$NV, DP = muts$DP, peaks = peaks) 
    genotypes = get_models(karyotype_combination = paste0(all_res_seg$k1[1],'-' ,all_res_seg$k2[1]),
                           purity = x$purity, ccf = all_res_seg$CCF_sim[1]) %>% filter(model_id == all_res_seg$model_id[1]) %>% 
      select(genotype_1,genotype_2) 
    g1 = genotypes %>% filter(!is.na(genotype_1)) %>% pull(genotype_1) %>% unique()
    g2 = genotypes %>% filter(!is.na(genotype_2)) %>% pull(genotype_2) %>% unique()
    b <- sum(log(BAF_LL(baf_obs = snps$BAF, n = snps$N.BAF,
                k1 = all_res_seg$k1[1], k2 = all_res_seg$k2[1],
                purity = patched_x$purity,
                ccf = all_res_seg$CCF_sim[1],
                g1 = g1, 
                g2 = g2)))
    
    d <- sum(log(DR_LL(dpr_obs = snps$DR,
               n = snps$N.BAF,
               k1 = all_res_seg$k1[1],
               k2 = all_res_seg$k1[1],
               purity = patched_x$purity, 
               ccf = all_res_seg$CCF_sim[1],
               ploidy=patched_x$ploidy)))
  }
}
  
  return(data.frame('vaf_ll'=v, 'baf_ll'=b, 'dr_ll'=d))
  
  
}

# patched_x : analysed cnaqc obj
# all_res_seg : results for a single segment
plot_simulation = function(patched_x, all_res_seg){
  
  n_snps = patched_x$snps %>% filter(segment_id==all_res_seg$segment_id[1]) %>% pull(N.BAF) %>% sum()
  n_points = patched_x$snps %>% filter(segment_id==all_res_seg$segment_id[1]) %>% pull(N.BAF) %>% length()
  n_points_vaf = patched_x$mutations  %>% filter(segment_id==all_res_seg$segment_id[1]) %>% pull(chr) %>% length()
  coverage= mean(patched_x$mutations$DP)
  st = 'aaaacc
        bbbbcc'
  
  muts = patched_x$mutations %>% filter(segment_id==all_res_seg$segment_id[1])
  snps = patched_x$snps %>% filter(segment_id==all_res_seg$segment_id[1])
  
  # Given solution name
  given_solution_name = paste0('Clonal ', all_res_seg$k1[1], ' purity ', patched_x$purity)
  # Correct solution name
  if (all_res_seg$k2[1] == ''){
    correct_solution_name = paste0('Clonal ', all_res_seg$k1[1], ' purity ', patched_x$purity)
  }else{
    correct_solution_name = paste0('Subclonal ', all_res_seg$k1[1],',', all_res_seg$k2[1], ' CCF ', all_res_seg$CCF_sim[1])
  }
  # Inferred solution name
  if (all_res_seg$k2_best[1] == ''){
    inferred_solution_name = paste0('Clonal ', all_res_seg$k1_best[1], ' purity ', patched_x$purity)
  }else{
    inferred_solution_name = paste0('Subclonal ', all_res_seg$k1_best[1],',', all_res_seg$k2_best[1], ' CCF ', all_res_seg$CCF_best[1])
  }
  
  
  given_baf = plot_baf_exp_vs_obs(patched_x, all_res_seg, n_snps, n_points, which='given') + ggtitle(given_solution_name)
  given_dr = plot_dr_exp_vs_obs(patched_x, all_res_seg, n_snps, n_points, which='given')
  given_vaf = plot_vaf_exp_vs_obs(patched_x, all_res_seg, n_points=n_points_vaf, coverage=coverage, which='given')
  given_solution_plot = patchwork::wrap_plots(given_baf, given_dr, given_vaf, design = st)
  # Given solution LL
  given_ll = get_likelihood_simulated_model(patched_x, all_res_seg, muts, snps, which='given') %>% mutate(id='given')
  
  # Correct solution plots
  correct_baf = plot_baf_exp_vs_obs(patched_x, all_res_seg, n_snps, n_points, which='correct') + ggtitle(correct_solution_name)
  correct_dr = plot_dr_exp_vs_obs(patched_x, all_res_seg, n_snps, n_points, which='correct')
  correct_vaf = plot_vaf_exp_vs_obs(patched_x, all_res_seg, n_points=n_points_vaf, coverage=coverage, which='correct')
  correct_solution_plot = patchwork::wrap_plots(correct_baf, correct_dr, correct_vaf, design = st)
  # Correct solution LL
  correct_ll = get_likelihood_simulated_model(patched_x, all_res_seg, muts, snps, which='correct') %>% mutate(id='correct')
  
  # Inferred solution
  inferred_baf = plot_baf_exp_vs_obs(patched_x, all_res_seg, n_snps, n_points, which='inferred') + ggtitle(inferred_solution_name)
  inferred_dr = plot_dr_exp_vs_obs(patched_x, all_res_seg, n_snps, n_points, which='inferred')
  inferred_vaf = plot_vaf_exp_vs_obs(patched_x, all_res_seg, n_points=n_points_vaf, coverage=coverage, which='inferred')
  inferred_solution_plot = patchwork::wrap_plots(inferred_baf, inferred_dr, inferred_vaf, design = st)
  # Inferred likelohood 
  inferred_ll = all_res_seg %>% select(vaf_ll,baf_ll,dr_ll) %>% mutate(id = 'inferred') 
  inferred_ll = inferred_ll[1,]
  
  ll_df = rbind(given_ll, correct_ll, inferred_ll) %>% pivot_longer(1:3) #reshape2::melt()
  
  ll_plot = ggplot(ll_df) +
    geom_col(aes(x= name, y=value, fill= id), position = "dodge")+theme_bw()+
    labs(x='',y='likelihood')+
    scale_fill_manual(values = c('correct'='forestgreen', 'given'='steelblue', 'inferred'='indianred'))
  
  compasison_plot = ggarrange(plotlist= list(given_solution_plot, correct_solution_plot, inferred_solution_plot), ncol= 1)
  
  final_plot = ggarrange(plotlist= list(compasison_plot, ll_plot), ncol= 2, heights = c(1, .5), widths = c(1, .5))
  
  final_plot
}



