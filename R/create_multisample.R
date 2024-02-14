#' Create multisample segmentation
#' 
#' @description Creates a m_CNAqc object starting from a named list of CNAqc objects on which the quality control has already been 
#' performed individually (ie: multiple samples coming from the same patient). Breakpoints from the original CNAqc objects are
#' used to define new breakpoints (and then segments) common across all the samples; in the resulting object, the mutations are 
#' remapped on the newly defined segments. 
#' The new m_CNAqc object can be used to perform new a quality control analysis with the different segmentation. 
#' 
#' @param cnaqc_objs named list of CNAqc objects, one per sample (each element must include at 
#' least "mutations", "cna", "reference" and "purity"). Names of the list must correspond to sample_id. 
#' 
#' @param cna_type can be "clonal", "subclonal" or both. Specifies which mutations must be considered in the creation of the 
#' m_CNAqc object. At the moment, the function is implemented for clonal mutations only.
#' 
#' @return A list of CNAqc objects, one per sample.
#' - `mutations` = mutations mapped on the new defined segments in the CNAqc format
#' - `cna` = table with the common segments in the CNAqc format
#' 
#' @examples # get better dataset! --> still to check
#' 
#' 
#' CNAqc_samples = lapply(names(sample_mutations), function(x) {
#'   CNAqc::init(mutations = snvs[[x]], 
#'            cna = cna[[x]], 
#'            purity = purity[[x]], 
#'            sample = x,
#'            ref = "GRCh38")
#'   })
#'
#' names(CNAqc_samples) = sapply(CNAqc_samples, function(x) {x$sample})
#'
#' multisamples = multisample_init(cnaqc_objs = CNAqc_samples)
#' 
#' multisamples



#####################

# initialize the multisample object with the new segmentation 

# multisample_init takes as input a list of cnaqc objects on which the quality control has already been done
# default type of mutations is clonal 

multisample_init <- function(cnaqc_objs, cna_type = "clonal") {
  
  #check if you are passing a list of cnaqc objs
  lapply(cnaqc_objs, function(x) {
    stopifnot(inherits(x, "cnaqc"))
  })
  
  # create the m_CNAqc object 
  multi_input = prepare_input_data_multiple(cnaqc_objs, cna_type)
  
  # define the output as a m_cnaqc object
  class(multi_input) <- "m_cnaqc"
  
  return(multi_input) 
  # the output is a m_cnaqc class object, in which each element of the list is a cnaqc object (one per sample) with all the classical attributes, but on which it has been performed the joint segmentation and mutations have been therefore remapped
}  
  

#####################################################################################
# internal functions

# prepare the data and set them in the correct format

prepare_input_data_multiple <- function(cnaqc_objs, cna_type) {
  
  #check if you are passing a list of cnaqc objs
  lapply(cnaqc_objs, function(x) {
    stopifnot(inherits(x, "cnaqc"))
  })
  
  # retrive information on breakpoints and define new segments (see the function for better explanation)
  # for the desidered type of mutations
  multi_cna <- join_segments(cnaqc_objs = cnaqc_objs, cna_type)
  
  # map the original mutations on the new segments for each segment
  multi_mutations <- lapply(names(multi_cna), function(x) {
    set_elements(cnaqc_obj = cnaqc_objs[[x]], new_cna = multi_cna[[x]], cna_type = cna_type)
    })
  
  names(multi_mutations) = names(multi_cna)
  
  # create a list with new cnaqc objects with the new segmentation and the mutations mapped on them
  multi_input = lapply(names(multi_cna), function(x) {
    init(mutations = multi_mutations[[x]], 
      cna = multi_cna[[x]], 
      purity = cnaqc_objs[[x]]$purity, 
      sample = x)
    })
  
  names(multi_input) <- lapply(multi_input, function(x) {x$sample}) %>% unlist()
  
  # assign m_cnaqc class to the output
  class(multi_input) <- "m_cnaqc"
  
  return(multi_input)

}
  
  
################################################################################################

# identification of new breakpoints
## the algorithm: 
## - binds together by rows all the cna elements of each sample included in the cnaqc objects list
## - iterates by chromosome, selecting the new breakpoints one chromosome at time 
## - defines the breakpoints by gathering together all the unique values included in the "from" and "to"
##    columns of the cna tibble, and then splitting it in two vectors of the same length, "new_from" and "new_to"
## - for each sample (considering one chr at time) searches if in the new range defined by "new_from" and "new_to" 
##    it was already defined a segment. If yes, it retrives information about it (ie: Major, minor, covRatio, BAF, etc),
##    changing the breakpoints, else fills all the columns (except for from, to and segment_id) with NA
## - pastes together everything in a unique tibble
## - removes segments that are not actually shared across samples (ie: if one sample does not have any cna event on a new segment the segment is eliminated)
##    --> !!! mutations that are mapped on not shared segments will not be included in the final object
## - divides the tibble by sample, returning a list of tibbles with the new defined segments, one per sample


########### getter of segment
## data = a tibble or a dataframe including all the original breakpoints from all the samples. 
## chr = chromosome on which iterate 
## sample = sample id
## new_from, new_to = newly defined breakpoints, used to retrive information about Major and minor allele copy number the
##                    genomic region defined by them
## keep_columns = names of columns that must be keeped in the new cna table. defined by removing a specific set of names  
##                from the col.names of the original tibble.

get_segment_info = function(data, chr, sample, new_from, new_to, keep_columns){
  data %>%
    dplyr::filter(sample_id == sample,
                  chr == !!chr,
                  from <= new_from,
                  to >= new_to) %>%
    select(all_of(keep_columns))
}

# segment definition

join_segments = function(cnaqc_objs, cna_type){
  
  # Row binded segments table (with sample specification)
  x = lapply(cnaqc_objs %>% names(), function(x){
    
    CNA(cnaqc_objs[[x]], type = cna_type) %>% 
      dplyr::mutate(sample_id = x)
    
  }) %>%
    do.call(bind_rows, .) %>%
    dplyr::select(sample_id, dplyr::everything())
  
  out = lapply(x$chr %>% unique(), function(chr) {
        
    # Chromosome-specific new breakpoints
    new_breakpoints = c(
      x %>%
        dplyr::filter(chr == !!chr) %>%
        dplyr::pull(from),
      x %>%
        dplyr::filter(chr == !!chr) %>%
        dplyr::pull(to)) %>% 
      unique() %>%
      sort()
    
    #  Separate new breakpoints into segment from and to values
    new_from = new_breakpoints[ !new_breakpoints == dplyr::last(new_breakpoints)] # last element will not be included in the from column 
    new_to = new_breakpoints[!new_breakpoints == dplyr::first(new_breakpoints)] # first element will not be included in the to column 
    
    
    # iterate on the new breakpoints to subset the cna piled up 
    lapply(new_breakpoints[-1] %>% seq_along(), function(i) {
      
      # iterate for each sample 
      lapply(x$sample_id %>% unique(), function(s) {
        
        # define which columns must be kept in the new table
        not_wanted = c("from", "to", "length", "size", "segment_id", "chr", "sample_id", "n")
        wanted = setdiff(colnames(x), not_wanted)
        
        
        # get the information for the sample in the new segment 
        tmp = get_segment_info(x, 
                               chr = chr, 
                               sample = s, 
                               new_from = new_from[i], 
                               new_to = new_to[i], 
                               keep_columns = wanted)
        
        # do some checking on the result
        if (nrow(tmp) == 0) { # there is no information on copy number on the new segment: insert NA as value of all the columns, except from, to, segment_id and sample_id
          
          tmp_v2 = rep(NA, ncol(tmp))  
          names(tmp_v2) = colnames(tmp)
          
          tmp = tmp_v2 %>% 
            tibble::as_tibble_row()
        }
        
        # create a tibble with the information on the new breakpoints and include the previously retrieved information 
        tidyr::tibble(chr = chr, 
               from = new_from[i], 
               to = new_to[i], 
               sample_id = s) %>% 
          dplyr::bind_cols(tmp) %>% 
          dplyr::mutate(segment_id = paste(chr, from, to, sep = ":")) 
        
      }) %>% do.call(bind_rows, .)
    }) %>% do.call(bind_rows, .)
  })  %>% do.call(bind_rows, .) # create a unique big tibble with the new segmentation
  
  
# remove all the segments that are not correctly shared across samples
  remove_segments = out %>% 
    dplyr::filter(is.na(Major)) %>% 
    dplyr::filter(is.na(minor)) %>% 
    dplyr::pull(segment_id) %>% 
    unique()
  
  all_segments = out %>% 
    pull(segment_id) %>% 
    unique()
  
  keep_segments = setdiff(all_segments, remove_segments)
  
  out = out %>% 
    dplyr::filter(segment_id %in% keep_segments)

  # split the cna table by sample id
  out_by_sample = lapply(out$sample_id %>% unique(), function(x) {
    out %>% 
      dplyr::filter(sample_id == x)
  })
  
  names(out_by_sample) = lapply(out_by_sample, function(x) {x$sample_id %>% unique}) %>% unlist()
    
  return(out_by_sample) 
  
  # returns a list with the new segmentation for each sample --> new cna of the m_cnaqc object
  
}

# map mutations on new segments
## uses CNAqc function "prepare_input_data" to map mutations on the newly defined segments.

set_elements <- function(cnaqc_obj, new_cna, cna_type) {

  # might need to define some errors
  
  initial_mutations = Mutations(cnaqc_obj, cna = cna_type)
  
  new_element = CNAqc:::prepare_input_data(initial_mutations, new_cna, cnaqc_obj$purity)

  remapped_mut = new_element$mutations

  return(remapped_mut)
}