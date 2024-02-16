#' printing errors, warnings and processing for multi_CNAqc objects and creation
#' 
#' 
#' 
#' 

# check if the input of multisample_init is actually a list of CNAqc objects and print the errors 
checking_input <- function(cnaqc_objs) {
  
  if (class(cnaqc_objs) != "list") {
    wrong_class_all = class(cnaqc_objs)
    
    cli::cli_abort(
      c("cnaqc_objs must be a list of CNAqc objects",
        "x" = "{.var cnaqc_objs} is not a list, you supplied a {.cls {class(cnaqc_objs)}} instead")
    )
  }
  
  if (unique(lapply(cnaqc_objs, class)) != "cnaqc") {
    wrong_class_elements = unique(unique(lapply(cnaqc_objs, class)))
    
    cli::cli_abort(
      c("{.var cnaqc_objs} must be a list of CNAqc objects",
        "x" = "Elements of {.var cnaqc_objs} are of class {.val {wrong_class_elements}} instead")
    )
  }
  
  if (names(cnaqc_objs) %>% is.null()) {
    
    cli::cli_abort(c("{.var cnaqc_objs} is not a named list", 
                   "x" = "Names of {.var cnaqc_objs} must correspond to the sample_id"))
  }
  
  if (lapply(cnaqc_objs, function(x) {x$sample}) %>% unique()%>% unlist %>% duplicated() %>% any() == TRUE) {
    
    cli::cli_abort(c("{.field sample} is the same in all CNAqc objects",
      "x"="Must provide a unique sample_id as {.field sample} in each CNAqc object"))
  }
  
  if ((lapply(cnaqc_objs, function(x) {x$peaks_analysis %>% is.null()}) %>% unlist %>% unique) %>% any() == TRUE) {
    
    cli::cli_abort(c("{.field peaks_analysis} is not present in all CNAqc objects",
                     "x" = "You should create a {.cls multi_CNAqc} object only after having performed the peak analysis on all samples"))
    
  } 
  
}
