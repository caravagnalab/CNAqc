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
#' m_CNAqc object. At the moment, the function is implemented for clonal mutations only. Default = "clonal"
#' 
#' @param QC_filter logical. Indicates wheter to filter or not for QC-passing mutations
#' 
#' @return a multi_CNAqc object. Every element of the object correspond to one of the sample, structured as follows: 
#' - `mutations_on_shared` = CNAqc object for the considered sample, containing all the information (mutations, cna, purity, etc) for
#'    mutations mapped on segments shared across all the samples;
#' - `mutations_on_private` = table including all the mutations, mapped on new segments that are not shared across all samples;
#' - `original_additional_info` = additional information (i.e.: CCF, peaks analysis, etc) from the original CNAqc object stored as a list (refers to the
#'    original segmentation!).
#' 
#' @export
#' 
#' @examples # get better dataset! --> still to check
#' 
#' 
#' 
#' CNAqc_samples = lapply(names(sample_mutations), function(x) {
#'   CNAqc::init(mutations = sample_mutations[[x]], 
#'            cna = cna[[x]], 
#'            purity = purity[[x]], 
#'            sample = x,
#'            ref = "GRCh38")
#'   })
#'
#' names(CNAqc_samples) = sapply(CNAqc_samples, function(x) {x$sample})
#' 
#' # perform peak analysis
#' 
#' CNAqc_samples = lapply(CNAqc_samples, function(x) {
#'   
#'   CNAqc::analyze_peaks(x, matching_strategy = 'closest')
#'   
#' })
#' 
#' create the multi_CNAqc object
#' 
#' multisamples = CNAqc::multisample_init(cnaqc_objs = CNAqc_samples)
#' 
#' multisamples
#' 
#' 
#' 
#' 
###################################################################################################################################################

# initialize the multisample object with the new segmentation 

# multisample_init takes as input a list of cnaqc objects on which the quality control has already been done
# default type of mutations is clonal 

multisample_init <- function(cnaqc_objs, 
                             cna_type = "clonal", 
                             QC_filter = TRUE) {
  
  cli::cli_h1("multi_CNAqc - Defining common segments")
  cat('\n')
  
  # perform some checkings on the input
  # - it is a list
  # - it contains all CNAqc objects
  # - it is a named list
  # - CNAqc objects all have the 'sample' field
  
  checking_input(cnaqc_objs)
  
  len = length(cnaqc_objs)
  cli::cli_alert_info("Selected CNA type: {.field {cna_type}}")
  cli::cli_alert_info("Found {.field {len}} CNAqc objects:")
  cli::cli_ul(names(cnaqc_objs)) 
  #cat('\n')
  
  # create the m_CNAqc object 
  
  cli::cli_h2("Building a {.cls multi_CNAqc} object")
  cat("\n")
  
  # retrive information on breakpoints and define new segments (see the function for better explanation)
  # for the desidered type of mutations
  
  cli::cli_rule("Selecting new segments")
  #cat("\n")
  
  multi_cna <- join_segments(cnaqc_objs = cnaqc_objs, cna_type, QC_filter)
  
  # map the original mutations on the new segments for each segment
  
  cli::cli_rule("Mapping mutations on new segments")
  cat("\n")
  
  multi_mutations <- lapply(names(multi_cna), function(x) {
    set_elements(cnaqc_obj = cnaqc_objs[[x]], new_cna_list = multi_cna[[x]], cna_type = cna_type, QC_filter = QC_filter)
  })
  
  names(multi_mutations) = names(multi_cna)
  
  # multi_input = prepare_input_data_multiple(cnaqc_objs, cna_type)
  
  # create a list with new cnaqc objects with the new segmentation and the mutations mapped on them
  
  cli::cli_h1("Defining the {.cls multi_CNAqc} object")
  
  # multi_input <- list()
  
  multi_input = lapply(names(multi_cna), function(x) {
    # list(
    shared = init(
      mutations = multi_mutations[[x]]$shared,
      cna = multi_cna[[x]]$shared,
      purity = cnaqc_objs[[x]]$purity,
      sample = x
    )#,
    private = multi_mutations[[x]]$private#,
    # original_additional_info = cnaqc_objs[[x]]$peaks_analysis
    # )
    
    original = setdiff(names(cnaqc_objs[[x]]), names(shared))
    other_info = lapply(original, function(o) {
      o = cnaqc_objs[[x]][[o]]
      })
    names(other_info) = original
    
    list(mutations_on_shared = shared, 
         mutations_on_private = private, 
         original_additional_info = other_info)
  })

  
  names(multi_input) <- lapply(multi_input, function(x) {x$mutations_on_shared$sample}) %>% unlist()
  class(multi_input) <- "m_cnaqc"
  # define the output as a m_cnaqc object
  
  
  if(exists("multi_input", inherits = F) & length(multi_input) == length(cnaqc_objs)) {
    cli::cli_h1("Ended")
    cli::cli_alert_success("{.cls multi_CNAqc} object created including all samples")
  }
  
  return(multi_input) 
  # the output is a m_cnaqc class object, in which each element of the list is a cnaqc object (one per sample) with all the classical attributes, but on which it has been performed the joint segmentation and mutations have been therefore remapped
}  
  

#####################################################################################
# internal functions

# prepare the data and set them in the correct format

# prepare_input_data_multiple <- function(cnaqc_objs, cna_type) {
# 
#   # retrive information on breakpoints and define new segments (see the function for better explanation)
#   # for the desidered type of mutations
#   multi_cna <- join_segments(cnaqc_objs = cnaqc_objs, cna_type)
#   
#   # map the original mutations on the new segments for each segment
#   multi_mutations <- lapply(names(multi_cna), function(x) {
#     set_elements(cnaqc_obj = cnaqc_objs[[x]], new_cna_list = multi_cna[[x]], cna_type = cna_type)
#     })
#   
#   names(multi_mutations) = names(multi_cna)
#   
#   # create a list with new cnaqc objects with the new segmentation and the mutations mapped on them
#   multi_input = lapply(names(multi_cna), function(x) {
#     init(mutations = multi_mutations[[x]], 
#       cna = multi_cna[[x]], 
#       purity = cnaqc_objs[[x]]$purity, 
#       sample = x)
#     })
#   
#   names(multi_input) <- lapply(multi_input, function(x) {x$sample}) %>% unlist()
#   
#   # assign m_cnaqc class to the output
#   class(multi_input) <- "m_cnaqc"
#   
#   return(multi_input)
# 
# }
  
  
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

join_segments = function(cnaqc_objs, cna_type, QC_filter){

  # Row binded segments table (with sample specification)
  x = lapply(cnaqc_objs %>% names(), function(x){
    
    CNA(cnaqc_objs[[x]], type = cna_type) %>% 
      dplyr::mutate(sample_id = x)
    
  }) %>%
    do.call(bind_rows, .) %>%
    dplyr::select(sample_id, dplyr::everything())
  
  if(QC_filter == TRUE) {
    x = x %>% 
      filter(QC_PASS == TRUE)
  } 
  
  out = lapply(x$chr %>% unique(), function(chr) {
    
    # cli::cli_alert_info("Iterating on {.field {chr}}")
    cli::cli_h2("{symbol$info} Iterating on {.field {chr}}")
    
    old_segments = sapply(x$sample_id %>% unique(), function(s) {
      x %>%
        dplyr::filter(sample_id == s) %>%
        dplyr::filter(chr == !!chr) %>%
        dplyr::pull(segment_id) %>%
        unique() %>%
        length()
    })
    
    cli::cli_alert_info("Number of original segments in individual {.cls CNAqc} objects in {.field {chr}}:")
    cli::cli_ul(paste(names(old_segments), old_segments, sep = " = "))
    cat("\n")
        
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
    
    cli::cli_alert_info("Found {.field {length(new_breakpoints)}} breakpoints")
    #cat("\n")
    
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
  
  # if(length(remove_segments) != 0) {
    
  cat("\n")
  cli::cli_alert_warning(
    "Found {.val {length(remove_segments)}} not shared segments"
  )
  
  all_segments = out %>%
    pull(segment_id) %>%
    unique()
  
  keep_segments = setdiff(all_segments, remove_segments)
  
  out_shared = out %>%
    dplyr::filter(segment_id %in% keep_segments)
  
  out_private = out %>%
    dplyr::filter(segment_id %in% remove_segments) %>%
    dplyr::filter(!is.na(minor) & !is.na(Major))
    
  cli::cli_alert_info("Shared segments across samples: {.val {out_shared %>% pull(segment_id) %>% unique() %>%  length()}}")
  
  # split the shared cna table by sample id
  out_shared_by_sample = lapply(out_shared$sample_id %>% unique(), function(x) {
    out_shared %>%
      dplyr::filter(sample_id == x)
  })
  
  # split the private cna table by sample id
  
  out_private_by_sample = lapply(out_private$sample_id %>% unique(), function(x) {
    if (nrow(out_private) == 0) {
      out_private %>%
        dplyr::add_row(chr = NA)
    } else {
      out_private %>%
        dplyr::filter(sample_id == x)
    }
  })
  
  if(length(out_private_by_sample) == 0) {
    lapply(x$sample_id %>% unique(), function(l) {
      out_private_by_sample[[l]] = rep(NA, ncol(out_private)) 
    })
  } else {
    names(out_private_by_sample) = lapply(out_private_by_sample, function(x) {x$sample_id %>% unique}) %>% unlist()
  }
  
  names(out_shared_by_sample) = lapply(out_shared_by_sample, function(x) {x$sample_id %>% unique}) %>% unlist()
  
  out_all_segments = lapply(x$sample_id %>% unique, function(x) {
    
    list(shared = out_shared_by_sample[[x]], private = out_private_by_sample[[x]])
    
  })
  
  names(out_all_segments) = lapply(out_all_segments, function(x) {x$shared$sample_id %>% unique()}) %>% unlist()
  
  cli::cli_alert_success("Obtained updated CNA table with shared and private segments per sample")  
  cat("\n")
  
  return(out_all_segments)
  
  # returns a list with the new segmentation for each sample --> new cna of the m_cnaqc object
  
}

# map mutations on new segments
## uses CNAqc function "prepare_input_data" to map mutations on the newly defined segments.

set_elements <- function(cnaqc_obj, new_cna_list, cna_type, QC_filter) {

  cli::cli_rule(
    crayon::bgCyan(crayon::white(cnaqc_obj$sample))
  )
  
  cat("\n")
  
  if(class(new_cna_list) != "list") {
    cli::cli_abort(c("Provided a {.cls {class(new_cna_list)}} object", 
                     "x" = "Input must be a list with the new segments"))
  }
  
  initial_mutations = Mutations(cnaqc_obj, cna = cna_type)
  
  if(QC_filter == TRUE) {
    initial_mutations = initial_mutations %>% 
      filter(QC_PASS == TRUE)
  }
  
  cli::cli_alert_info("Found {.val {nrow(initial_mutations)}} mutations in the original {.cls CNAqc} object")
  cat("\n")
  
  cli::cli_rule(paste(crayon::cyan("{symbol$info}"), "Mapping mutations on shared segments"))
  mutations_shared_segments = CNAqc:::prepare_input_data(mutations = initial_mutations, cna = new_cna_list$shared, tumour_purity = cnaqc_obj$purity)
  remapped_mut_shared = mutations_shared_segments$mutations
  
  cat("\n")
  cli::cli_rule(paste(crayon::cyan("{symbol$info}"), "Mapping mutations on private segments"))
  
  if(new_cna_list$private %>% is.na() %>% all() == TRUE) {
    
    cli::cli_alert_warning("No private segments found")
    remapped_mut = list(shared = remapped_mut_shared, private = NA)
    
  } else {
    
    mutations_private_segments = CNAqc:::prepare_input_data(mutations = initial_mutations, cna = new_cna_list$private, tumour_purity = cnaqc_obj$purity)
    remapped_mut_private = mutations_private_segments$mutations
    remapped_mut = list(shared = remapped_mut_shared, private = remapped_mut_private)
    cat("\n")
  }
  
  return(remapped_mut)
}