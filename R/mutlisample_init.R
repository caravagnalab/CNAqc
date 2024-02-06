#' multisample utilities
#' 
#' @description
#' In order to obtain a "joint table" including all the different samples of the same patient the
#' existing CNAqc objects (one per sample) must be merged, obtaining a new multi-sample CNAqc object
#' (mCNAqc) which includes all the information about the mutations and a new common segmentation
#' between all samples. 
#' 
#' @param cnaqc_objs named list of CNAqc objects, one per sample (each element must include at 
#' least "mutations", "cna", "reference" and "purity"). Names of the list must correspond to sample_id 
#' 
#' @examples # get better dataset!
#' x1 = init(
#'   mutations = example_dataset_CNAqc$mutations %>%  sample_n(1000),
#'   cna = example_dataset_CNAqc$cna %>% sample_n(100),
#'   purity = example_dataset_CNAqc$purity,
#'   sample = "xxx",
#'   ref = example_dataset_CNAqc$reference,
#' )

#' x2 = init(
#'   mutations = example_dataset_CNAqc$mutations %>%  sample_n(1000),
#'   cna = example_dataset_CNAqc$cna %>% sample_n(100),
#'   purity = example_dataset_CNAqc$purity,
#'   sample = "yyy",
#'   ref = example_dataset_CNAqc$reference,
#' )
#' 

#' cnaqc_objs = list(x1,x2)
#' names(cnaqc_objs) = c("s1", "s2")
#' 
#' mCNAqc = multisample_init(cnaqc_objs = cnaqc_objs)
#' 
#' @return 
#' A list of four elements: 
#' - `joint_segments` = table with the join segmentation. output of the join_segments function
#' - `joint_table` = table with mutations mapped on the new common segments 
#' - `metadata` = a tibble including additional information for each sample, such as reference genome and purity
#' - `original_data` = a list in which each element is the original CNAqc object for every sample

#####################

# initialize the multisample object

multisample_init <- function(cnaqc_objs) {
  
  #check if you are passing a list of cnaqc objs
  lapply(cnaqc_objs, function(x) {
    stopifnot(inherits(x, "cnaqc"))
  })
  
  output <- list()
  
  output[["joint_segments"]] <- join_segments(cnaqc_objs = cnaqc_objs)
  output[["joint_table"]] <-
    join_mutations(cnaqc_objs = cnaqc_objs, joint_segments = output[["joint_segments"]])
  output[["metadata"]] <- tibble(sample = cnaqc_objs %>% names(),
                                 purity = sapply(cnaqc_objs, function(x)
                                   x$purity))
  output[["original_data"]] <- cnaqc_objs
  
  return(output)
  
}

## internal functions

# perform joint segmentation

join_segments = function(cnaqc_objs){
  # Row binded segments table (with sample specification)
  x = lapply(cnaqc_objs %>% names(), function(x){
    cnaqc_objs[[x]]$cna$sample = x
    return(cnaqc_objs[[x]]$cna)
  }) %>%
    do.call(rbind, .) %>%
    dplyr::select(sample, dplyr::everything())
  
  out = lapply(x$chr %>% unique(), function(chr) {
    
    # Chromosome-specific new breakpoints
    new_breakpoints = c(
      x %>%
        dplyr::filter(chr == !!chr) %>%
        dplyr::pull(from),
      x %>%
        dplyr::filter(chr == !!chr) %>%
        dplyr::pull(to)
    ) %>% unique() %>%  sort()
    
    #  Separate new breakpoints into segment from and to values
    new_from = new_breakpoints[-length(new_breakpoints)]
    new_to = new_breakpoints[-1]
    
    lapply(new_from %>% seq_along(), function(i) {
      lapply(x$sample %>% unique(), function(s) {
        
        # Get sample-specific segment quantities (Major and minor copy numbers, covRatio)
        get_segment_info = function(x, what){
          x %>%
            dplyr::filter(sample == s,
                          chr == !!chr,
                          from <= new_from[i],
                          to >= new_to[i]) %>%
            dplyr::pull(what)
        }
        
        # keep some of the original columns
        
        not_wanted = c("from", "to", "length", "size", "segment_id", "chr", "sample")
        wanted = setdiff(colnames(x), not_wanted)
        
        tmp = sapply(wanted, function(w) get_segment_info(x, what = w))
        
        # check the structure of the result
        if (class(tmp) == list) {
          tmp = tmp %>% 
            as_tibble
        }
        
        if (nrow(tmp) == 0) {
          
        }
        if (nrow(tmp) == 1) {
          tmp = tmp %>% as_tibble_row
        } 
        else {
          tmp = tmp %>% as.tibble()
        }
        #  }
        
        tibble(chr = chr, 
               from = new_from[i], 
               to = new_to[i], 
               sample = s
               ) %>% bind_cols(tmp)
        
        #Major = get_segment_info(x, what = "Major")
        
        #minor = get_segment_info(x, what = "minor") 
        
        #  Return joint segments table
        tibble(
          chr = chr,
          from = new_from[i],
          to = new_to[i],
          sample = s,
          Major = ifelse(length(Major) == 0, NA, Major), #  Replace missing values with NA
          minor = ifelse(length(Major) == 0, NA, minor), #  Replace missing values with NA
          covRatio = ifelse(length(covRatio) == 0, NA, covRatio) #  Replace missing values with NA
        ) %>%
          mutate(segment_id = paste(chr, from, to, sep = ":"))
      }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .)
  })  %>% do.call(rbind, .)
  
  return(out)
}

# map the present mutations on the new segments

join_mutations = function(cnaqc_objs, joint_segments){
  
  # create the join table with all the mutations for each sample 
  x = lapply(cnaqc_objs %>% names(), function(x){
    cnaqc_objs[[x]]$mutations$sample = x
    return(cnaqc_objs[[x]]$mutations)
  }) %>%
    do.call(rbind, .) %>%
    dplyr::select(sample, dplyr::everything()) 
  
  # map the mutations on the new joint segments
  
  out = lapply(x$sample %>% unique(), function(s) {
    which_segments = joint_segments %>% dplyr::filter(sample == s) #get the new segments for each sample
    
    which_mutations = x %>%
      dplyr::filter(sample == s) %>%
      dplyr::mutate(segment_id = NA,
                    karyotype = NA) #initialize segment_id and karyotype columns (or override if they already exist)
    
    # iterate over the new segments and map the mutations on them 
    lapply(1:nrow(which_segments), function(i) {
      mappable = which(
        which_mutations$chr == which_segments$chr[i] &
          which_mutations$from >= which_segments$from[i] &
          which_mutations$to <= which_segments$to[i]
      )
      
      # fill karyotype and segment_id columns with the information for each mapped mutation. If the mutation does not map any
      # segment, it will be flagged as NA
      which_mutations$karyotype[mappable] = paste0(which_segments$Major[i], ':', which_segments$minor[i]) 
      which_mutations$segment_id[mappable] = which_segments$segment_id[i]
      
      return(which_mutations)
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .) %>%
    dplyr::mutate(karyotype = ifelse(karyotype == "NA:NA", NA, karyotype)) 
  
  return(out)
}
