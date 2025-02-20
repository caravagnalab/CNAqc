#' get sample names from CNAqc and mCNAqc objects
#'
#' @param x a CNAqc or a mCNAqc object.
#'
#' @return a vector with sample names of CNAqc objects. 
#' @export
#'
#' @examples
#' \dontrun{get_sample_name(x)}
#'
get_sample_name <- function(x) {
  
  if (class(x) == "m_cnaqc") {
    
    lapply(x$cnaqc_obj_new_segmentation, function(y) {
      
      y$sample
      
    }) %>% unlist() %>% unname()
    
  } else if (class(x) == "cnaqc") {
    
    x$sample
    
  } else {
   
    wrong_class_all = class(x)
    
    cli::cli_abort(
      c("must provide a {.field m_cnaqc} object",
        "x" = "{.var x} is a {.cls {class(x)}}")
    )
  }
}

#############################################
#' get desidered samples
#'
#' @param m_cnaqc_obj mCNAqc object. 
#' @param sample list of one or more samples, must correspond to the names of each element in the mCNAqc object
#' @param which_obj one of \code{shared} or \code{original}.
#'
#' @return a list of CNAqc objects corresponding to the desidered samples.
#' @export
#'
#' @examples 
#' \dontrun{get_sample(x, sample = get_sample_names(x), which_obj = "shared") }
#' 
get_sample <- function(m_cnaqc_obj,
                       sample,
                       which_obj) {
  if (class(m_cnaqc_obj) != "m_cnaqc") {
    wrong_class_all = class(m_cnaqc_obj)
    
    cli::cli_abort(
      c("cnaqc_objs must be a {.field m_cnaqc} object",
        "x" = "{.var m_cnaqc_obj} is a {.cls {class(m_cnaqc_obj)}}")
    )
  }
  
  consented_obj <- c("shared", "original")
  
  if ((which_obj %in% consented_obj) == FALSE) {
    
    cli::cli_abort("{.var which_obj} must be one of {.val shared} or {.val original}")
  
  }
  
  # define the element names
  
  if (which_obj == "original") {
    
    type = "original_cnaqc_objc"
    
    # check if the original cnaqc obj exist
    
    check_or = any(names(m_cnaqc_obj) == type)

    if (check_or == FALSE) {
      cli::cli_abort(c("mCNAqc object was build without keeping original CNAqc objects"), 
                     "x" = "It is not possible to retrieve the required samples")
    } else {
      
      cli::cli_h1("Retrieving original {.cls CNAqc} objects")
      cnaqc_samples = m_cnaqc_obj[[type]][sample]
      
    }
    
  } else {
    
    type = "cnaqc_obj_new_segmentation"
    
    cli::cli_h1("Retrieving {.cls CNAqc} objects with the new segmentation")
    
    cnaqc_samples = m_cnaqc_obj[[type]][sample]
  }
  
  return(cnaqc_samples)
}

#' Get mCNAqc statistics
#'
#' @param m_cnaqc_obj 
#'
#' @return a tibble
#' @export
#'
#' @examples
#' \dontrun{get_mCNAqc_stats(x)}
#' 
get_mCNAqc_stats <- function(m_cnaqc_obj){
  stats = m_cnaqc_obj[["m_cnaqc_stats"]]
  return(stats)
}