#' Extract mutations from common and private segments from multi_CNAqc object. 
#' 
#' @description 
#' This function extracts the all the mutations in CNAqc format of desired samples from a multi_CNAqc object, 
#' merging together those that belong to common segments and those private to the sample.
#' 
#' @param m_cnaqc_obj multi_CNAqc object (must have previously run \code{"multisample_init"})
#' @param sample desired samples for which you want to extract the mutations. Default \code{"all"}, extracts mutations from all samples
#' in the multi_CNAqc object. 
#' @param cna \code{"clonal"} for clonal CNAs, \code{"subclonal"} for subclonal CNAs.
#' @param type \code{"SNV"} for single-nucleotide variants, \code{"indel"} for insertion-deletions.
#'
#' @return a named list of tibbles with all the mutations.
#' @export
#'
#' @examples
#'  
#' all_segments_mut = All_segments_mutations(example_multisample, 
#'                                 sample = "all", 
#'                                 cna = c("clonal", "subclonal"), 
#'                                 type = c("SNV", "indel"))
#' print(all_segments_mut)

All_segments_mutations = function(m_cnaqc_obj, 
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

#' Extract mutations in common segments from multi_CNAqc object. 
#' 
#' @description 
#' This function extracts the all the mutations belonging to common segments in CNAqc format of 
#' desired samples from a multi_CNAqc object.
#' 
#' @param m_cnaqc_obj multi_CNAqc object (must have previously run \code{"multisample_init"})
#' @param sample desired samples for which you want to extract the mutations. Default \code{"all"}, extracts mutations from all samples
#' in the multi_CNAqc object. 
#' @param cna \code{"clonal"} for clonal CNAs, \code{"subclonal"} for subclonal CNAs.
#' @param type \code{"SNV"} for single-nucleotide variants, \code{"indel"} for insertion-deletions.
#'
#' @return a named list of tibbles with all the mutations.
#' @export
#'
#' @examples
#'  
#' shared_segments_mut = Shared_segments_mutations(example_multisample, 
#'                                 sample = "all", 
#'                                 cna = c("clonal", "subclonal"), 
#'                                 type = c("SNV", "indel"))
#' print(shared_segments_mut)

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