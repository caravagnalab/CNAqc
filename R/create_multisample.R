#' Create multisample segmentation
#' 
#' @description Creates a mCNAqc object starting from a named list of CNAqc objects on which the quality control has already been 
#' performed individually (ie: multiple samples coming from the same patient). Breakpoints from the original CNAqc objects are
#' used to define new breakpoints (and then segments) common across all the samples; in the resulting object, the mutations are 
#' remapped on the newly defined segments. 
#' The new m_CNAqc object can be used to perform new a quality control analysis with the different segmentation. 
#' 
#' @param cnaqc_objs named list of CNAqc objects, one per sample (each element must include at 
#' least "mutations", "cna", "reference" and "purity"). Names of the list must correspond to sample_id. 
#' 
#' @param QC_filter logical. Indicates whether to filter or not for QC-passing mutations (NB: must have performed the peak analysis before). 
#' Default is \code{TRUE}.
#' 
#' @param keep_original logical. Indicates whether to keep or not original CNAqc objects, default is \code{TRUE}
#' 
#' @param discard_private logical. Indicates whether to keep or not mutations falling on private segments, default is \code{FALSE}.
#' 
#' @return a mCNAqc object, structured as follows: 
#' - `cnaqc_obj_new_segmentation` = list of CNAqc objects for all the samples created using the new segmentation;
#' - `original_cnaqc_objc` = original CNAqc object, only if the argument \code{keep_original} is set to TRUE;
#' - `original_additional_info` = additional information from the original CNAqc objects, such peaks analysis results of CCF estimation. Only if the argument \code{keep_original} is set to FALSE;
#' - `m_cnaqc_stats` = information related to number of mutations and segments before and after the creationn of the mCNAqc object. 
#' 
#' @export
#' 
#' @examples 
#' \dontrun{
#' multisamples = CNAqc::multisample_init(cnaqc_objs = CNAqc_samples, 
#'                  QC_filter = TRUE, 
#'                  keep_original = FALSE, 
#'                  discard_private = FALSE)
#' 
#' multisamples
#' }
#' 
###################################################################################################################################################

# initialize the multisample object with the new segmentation 

# multisample_init takes as input a list of cnaqc objects on which the quality control has already been done
# default type of mutations is clonal 

multisample_init <- function(cnaqc_objs, 
                             # cna_type = "clonal",
                             QC_filter = TRUE, 
                             keep_original = TRUE, 
                             discard_private = FALSE) {
  
  cli::cli_h1("mCNAqcqc - Defining common segments")
  cat('\n')
  
  # perform some checks on the input
  # - it is a list
  # - it contains all CNAqc objects
  # - it is a named list
  # - CNAqc objects all have the 'sample' field
  
  checking_input(cnaqc_objs)
  
  len = length(cnaqc_objs)
  # cli::cli_alert_info("Selected CNA type: {.field {cna_type}}")
  cli::cli_alert_info("Found {.field {len}} CNAqc objects:")
  cli::cli_ul(names(cnaqc_objs)) 
  #cat('\n')
  
  # create the m_CNAqc object 
  
  cli::cli_h2("Building a {.cls mCNAqcqc} object")
  cat("\n")
  
  # retrive information on breakpoints and define new segments (see the function for better explanation)
  # for the desidered type of mutations
  
  cli::cli_rule("Selecting new segments")
  #cat("\n")
  
  multi_cna <- join_segments(cnaqc_objs = cnaqc_objs, 
                             # cna_type, 
                             QC_filter, 
                             keep_original)
  
  # map the original mutations on the new segments for each segment
  
  cli::cli_rule("Mapping mutations on new segments")
  cat("\n")
  
  multi_mutations <- lapply(names(multi_cna), function(x) {
      set_elements(cnaqc_obj = cnaqc_objs[[x]], 
                 new_cna_list = multi_cna[[x]], 
                 # cna_type = cna_type, 
                 QC_filter = QC_filter#, 
                 # keep_original
                 )
  })
  
  names(multi_mutations) = names(multi_cna)
  
  if(discard_private == TRUE) {
  
    # take only mutations on identical positions across all the samples
    cat("\n")
    cli::cli_rule("Collecting only mutations on shared positions")
    
    shared_mut = lapply(multi_mutations, function(x) {
      x %>%
        dplyr::mutate(pos = paste(chr, from, to, sep = ":"))
    }) %>% dplyr::bind_rows(.)
    
    n_samples = shared_mut$Indiv %>% unique() %>% length()
    tot_n_mut = nrow(shared_mut)
    shared_mut = shared_mut %>%
      dplyr::group_by(pos) %>%
      dplyr::filter(n() == n_samples) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(pos = NULL)
    
    final_n_mut = nrow(shared_mut)
    removed_n_mut = tot_n_mut - final_n_mut
    
    cli::cli_alert_info(
      c(
        "Found {.val {tot_n_mut}} mutations mapping common segments across {.val {n_samples}} samples. \n",
        "Removing {.val {removed_n_mut}} mutations on not shared positions, keeping {.val {final_n_mut}} mutations"
      )
    )
    cat("\n")
    
    multi_mutations = lapply(names(multi_mutations), function(x) {
      shared_mut %>%
        filter(Indiv == x)
    })
    names(multi_mutations) = names(multi_cna)
  }
  
  # create a list with new cnaqc objects with the new segmentation and the mutations mapped on them
  
  cli::cli_h1("Defining the {.cls mCNAqcqc} object")
  cat("\n")
  
  # creating the output 
  m_cnaqc_res <- list()
  
  # define the output as a m_cnaqc object
  class(m_cnaqc_res) <- "m_cnaqc"
  
  # creating the elements of the mCNAqc object
    # list of CNAqc obj with the new segments
  
  new_segmentation_cnaqc <- lapply(names(multi_cna), function(x) {
    #print(x)
    CNAqc::init(
      mutations = multi_mutations[[x]],
      cna = multi_cna[[x]]$shared,
      purity = cnaqc_objs[[x]]$purity,
      sample = x
    )})
  
  names(new_segmentation_cnaqc) <- sapply(new_segmentation_cnaqc, function(x) {x$sample})
  
  # multi_input = lapply(names(multi_cna), function(x) {
  #   #print(x)
  #   shared = CNAqc::init(
  #     mutations = multi_mutations[[x]],
  #     cna = multi_cna[[x]]$shared,
  #     purity = cnaqc_objs[[x]]$purity,
  #     sample = x
  #   )
  
  m_cnaqc_res$cnaqc_obj_new_segmentation <- new_segmentation_cnaqc
  
  if (keep_original == TRUE) {
    # private = multi_mutations[[x]]$private
    # list(mutations_on_shared = shared,
         # mutations_on_private = private,
         # original = cnaqc_objs[[x]]) %>% return()
    
    m_cnaqc_res$original_cnaqc_objc <- cnaqc_objs
    
  } else {
     
    original_elements <- lapply(cnaqc_objs, function(x) {names(x)}) %>% unlist() %>% unique()
    conserved_elements <- lapply(m_cnaqc_res$cnaqc_obj_new_segmentation, function(x) {names(x)}) %>% unlist() %>% unique()
    
    to_include = setdiff(original_elements, conserved_elements)
    
    other_info = lapply(cnaqc_objs, function(x) {x[to_include]})
    m_cnaqc_res$original_additional_info <- other_info
    
    # other_info = lapply(original, function(o) {
    #   cnaqc_objs[[x]][[o]]
    # })
    # names(other_info) = original
    # list(mutations_on_shared = shared,
    #      # mutations_on_private = private,
    #      original_additional_info = other_info) %>% return()
  }
  
  
  # names(multi_input) <- lapply(multi_input, function(x) {x$mutations_on_shared$sample}) %>% unlist()
  
  cli::cli_h2("Creating mCNAqc stats")
  
  n_new_mut <- sapply(m_cnaqc_res$cnaqc_obj_new_segmentation, function(x) {x$n_mutations})
  n_or_mut <- sapply(cnaqc_objs, function(x) {x$n_mutations})
  n_new_cna <- sapply(m_cnaqc_res$cnaqc_obj_new_segmentation, function(x) {x$n_cna})
  n_or_cna <- sapply(cnaqc_objs, function(x) {x$n_cna})
  
  stats_mcnaqc <- data.frame(
    n_mutations_original_segmentation = n_or_mut,
    n_mutations_new_segmentation = n_new_mut, 
    n_cna_original_segmentation = n_or_cna, 
    n_cna_new_segmentation = n_new_cna
  )
  
  m_cnaqc_res$m_cnaqc_stats <- stats_mcnaqc
  
  if(exists("m_cnaqc_res", inherits = F) & length(m_cnaqc_res$cnaqc_obj_new_segmentation) == length(cnaqc_objs)) {
    cli::cli_h1("Ended")
    cli::cli_alert_success("{.cls mCNAqc} object created including all samples")
  }
  
  return(m_cnaqc_res) 
  # the output is a m_cnaqc class object, in which each element of the list is a cnaqc object (one per sample) 
  # with all the classical attributes, but on which it has been performed the joint segmentation and mutations 
  # have been therefore remapped
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

join_segments = function(cnaqc_objs, 
                         # cna_type, 
                         QC_filter, 
                         keep_original){

  # Row binded segments table (with sample specification)
  x = lapply(cnaqc_objs %>% names(), function(x){
    
    CNA(cnaqc_objs[[x]]#, 
        # type = cna_type
        ) %>% 
      dplyr::mutate(sample_id = x)
    
  }) %>%
    do.call(bind_rows, .) %>%
    dplyr::select(sample_id, dplyr::everything())
  
  if(QC_filter == TRUE) {
    x = x %>% 
      filter(QC_PASS == TRUE)
  } 
  
  out = lapply(x$chr %>% unique(), function(chr) {
    
    cli::cli_alert_info("Iterating on {.field {chr}}")
    cli::cli_h2("Iterating on {.field {chr}}")
    
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
    
  all_segments = out %>%
    pull(segment_id) %>%
    unique()
  
  keep_segments = setdiff(all_segments, remove_segments)
  
  out_shared = out %>%
    dplyr::filter(segment_id %in% keep_segments)
    
  cli::cli_alert_info("Shared segments across samples: {.val {out_shared %>% pull(segment_id) %>% unique() %>%  length()}}")
  
  # split the shared cna table by sample id
  out_by_sample = lapply(out_shared$sample_id %>% unique(), function(x) {
    list(shared = out_shared %>%
      dplyr::filter(sample_id == x))
  })
  names(out_by_sample) = lapply(out_by_sample, function(x) {lapply(x, function(y) {y$sample_id %>% unique}) %>% unlist()}) %>% unlist() %>% unname()
  
  # select mutations on private segments
  
  if (keep_original == TRUE) {
    
    cli::cli_alert_warning("param {.arg keep_original} is set to {.var TRUE}. Original CNAqc object will be kept.")
    
    # out_private = out %>%
    #   dplyr::filter(segment_id %in% remove_segments) %>%
    #   dplyr::filter(!is.na(minor) & !is.na(Major))
    
    # split the private cna table by sample id

    # out_private_by_sample = lapply(out_private$sample_id %>% unique(), function(x) {
    #   if (nrow(out_private) == 0) {
    #     out_private %>%
    #       dplyr::add_row(chr = NA)
    #   } else {
    #     out_private %>%
    #       dplyr::filter(sample_id == x)
    #   }
    # })
    # 
    # if (length(out_private_by_sample) == 0) {
    #   lapply(x$sample_id %>% unique(), function(l) {
    #     out_private_by_sample[[l]] = rep(NA, ncol(out_private))
    #   })
    # } else {
    #   names(out_private_by_sample) = lapply(out_private_by_sample, function(x) {
    #     x$sample_id %>% unique
    #   }) %>% unlist()
    # }
    # 
    # out_by_sample = lapply(names(out_by_sample), function(x) {
    #   out_by_sample[[x]] = list(shared = out_by_sample[[x]]$shared, private = out_private_by_sample[[x]])
    # })
    
  } else {
    cat("\n")
    cli::cli_alert_warning(c(
      "param {.arg keep_original} is set to {.var FALSE}.",
      "Found {.val {length(remove_segments)}} not shared segments. Original CNAqc object will not be saved"
      )
    )
  }
  
  names(out_by_sample) = lapply(out_by_sample, function(x) {x$shared$sample_id %>% unique()}) %>% unlist()
  
  cli::cli_alert_success("Obtained updated CNA table with shared segments across samples")  
  cat("\n")
  
  return(out_by_sample)
  
  # returns a list with the new segmentation for each sample --> new cna of the m_cnaqc object
  
}

# map mutations on new segments
## uses CNAqc function "prepare_input_data" to map mutations on the newly defined segments.

set_elements <- function(cnaqc_obj, 
                         new_cna_list, 
                         # cna_type, 
                         QC_filter #, 
                         # keep_original
                         ) {

  cli::cli_rule(
    crayon::bgCyan(crayon::white(cnaqc_obj$sample))
  )
  
  cat("\n")
  
  if(class(new_cna_list) != "list") {
    cli::cli_abort(c("Provided a {.cls {class(new_cna_list)}} object", 
                     "x" = "Input must be a list with the new segments"))
  }
  
  initial_mutations = CNAqc::Mutations(cnaqc_obj
                                       # cna = cna_type
                                       )
  
  if(QC_filter == TRUE) {
    initial_mutations = initial_mutations %>% 
      filter(QC_PASS == TRUE)
  }
  
  cli::cli_alert_info("Found {.val {nrow(initial_mutations)}} mutations in the original {.cls CNAqc} object")
  cat("\n")
  
  cli::cli_rule("Mapping mutations on shared segments")
  mutations_shared_segments = CNAqc:::prepare_input_data(mutations = initial_mutations, cna = new_cna_list$shared, tumour_purity = cnaqc_obj$purity)
  remapped_mut = mutations_shared_segments$mutations
  
  # if (keep_original == TRUE) {
  #   cat("\n")
  #   cli::cli_rule(
  #     paste(
  #       crayon::cyan("{symbol$info}"),
  #       "Storing mutations on private segments in the {.cls mCNAqc} object"
  #     )
  #   )
  #   
  #   if (new_cna_list$private %>% is.na() %>% all() == TRUE) {
  #     #
  #     cli::cli_alert_warning("No private segments found")
  #     remapped_mut = list(shared = remapped_mut_shared, private = NA)
  #     #
  #   } else {
  #     #
  #     mutations_private_segments = CNAqc:::prepare_input_data(
  #       mutations = initial_mutations,
  #       cna = new_cna_list$private,
  #       tumour_purity = cnaqc_obj$purity
  #     )
  #     remapped_mut_private = mutations_private_segments$mutations
  #     remapped_mut = list(shared = remapped_mut_shared, private = remapped_mut_private)
  #     cat("\n")
  #   }
  # } else {
    # remapped_mut = list(shared = remapped_mut_shared)
  # }
  
  return(remapped_mut)
}