# defining getters for m_CNAqc object

# getter of all mutations of a sample (private + shared)

Private_segments_mutations = function(m_cnaqc_obj, 
                            sample = "all", 
                            cna = c("clonal", "subclonal"), 
                            type = c("SNV", "indel")) {
  
  if (class(m_cnaqc_obj) != "m_cnaqc") {
    wrong_class_all = class(m_cnaqc_obj)
    
    cli::cli_abort(
      c("cnaqc_objs must be a {.field m_cnaqc} object",
        "x" = "{.var m_cnaqc_obj} is a {.cls {class(m_cnaqc_obj)}}")
    )
  }
  
  # get the mutations for required samples
  
  if (sample == "all") {
    samples = names(m_cnaqc_obj)
  } else {
    samples = sample
  }
  
  sample_mutations = lapply(samples, function(x) {
    
    mut_on_private_segments = m_cnaqc_obj[[x]][["private_mutations"]]
    
    mut_on_shared_segments = Mutations(m_cnaqc_obj[[x]][["shared_mutations"]], cna = cna, type = type)
    
    all_mutations = dplyr::bind_rows(mut_on_private_segments, mut_on_shared_segments) #%>%
      # dplyr::arrange(desc(chr))
    return(all_mutations)
  })
  
  names(sample_mutations) = samples
  
  return(sample_mutations)
}

# getter of all common mutations across samples (shared a1 + shared a2)

Shared_segments_mutations = function(m_cnaqc_obj,
                                     sample = "all",
                                     cna = c("clonal", "subclonal"),
                                     type = c("SNV", "indel")) {
  if (class(m_cnaqc_obj) != "m_cnaqc") {
    wrong_class_all = class(m_cnaqc_obj)
    
    cli::cli_abort(
      c("cnaqc_objs must be a {.field m_cnaqc} object",
        "x" = "{.var m_cnaqc_obj} is a {.cls {class(m_cnaqc_obj)}}")
    )
  }
  
  # get the mutations for required samples
  
  if (sample == "all") {
    samples = names(m_cnaqc_obj)
  } else {
    samples = sample
  }
  
  mut_on_shared_segments = lapply(samples, function(x) {
    Mutations(m_cnaqc_obj[[x]][["shared_mutations"]], cna = cna, type = type)
  }) %>% 
    dplyr::bind_rows(.)
   
  return(mut_on_shared_segments)
  
}
