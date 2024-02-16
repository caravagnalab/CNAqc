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
  
  # if (lapply(cnaqc_objs, function(x) {exists(x$sample)})) add checking that the sample field exists in each cnaqc object
}


