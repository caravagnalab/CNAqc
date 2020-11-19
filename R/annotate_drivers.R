#' Title
#'
#' @param x 
#' @param drivers_list 
#' @param function_column 
#' @param gene_column 
#'
#' @return
#' @export
#'
#' @examples
annotate_drivers = function(x, 
                            drivers_list = CNAqc::DriverDBv3,
                            function_entries = c('exonic', "Exonic"),
                            function_column = 'ANNOVAR_FUNCTION',
                            gene_column = "GENE"
)
{
  stopifnot(inherits(x, 'cnaqc'))
  
  cli::cli_alert_info("Looking for exonic mutations using column {gene_column} and {function_column}: entries will be {function_entries}.")
  if(!(gene_column %in% colnames(x$snvs))) stop("Gene column ", gene_column, " is not available in the data:\n\t", paste(colnames(x$snvs), collapse = ', '), '.')
  if(!(function_column %in% colnames(x$snvs))) stop("Function column ", function_column, " is not available in the data:\n\t", paste(colnames(x$snvs), collapse = ', '), '.')
  
  to_find = x$snvs[[gene_column]][x$snvs[[function_column]] %in% function_entries] %>% sort %>% unique
  
  cli::cli_alert_info("Found exonic mutations in {length(to_find)} genes")
  print(to_find)
  
  found = drivers_list %>% 
    filter(driver_gene %in% to_find) %>% 
    pull(driver_gene) %>% 
    unique
  
  cli::cli_alert_info("{length(found)} found the input list of {length(drivers_list$driver_gene) %>% unique()} genes.")
  if(length(found) > 0) print(found)
  
  if(found == 0) return(x)
  
  which_driver = intersect(
    which(x$snvs[[gene_column]] %in% drivers_list$driver_gene),
    which(x$snvs[[function_column]] %in% function_entries)
  )
  
  x$snvs$is_driver = FALSE
  x$snvs$is_driver[which_driver] = TRUE
  
  x$snvs %>% 
    filter(is_driver) %>% 
    print
  
  cli::cli_alert_info("Reassembling a new CNAqc object with these drivers")
  x = CNAqc::init(
    x$snvs,
    cna = x$cna,
    purity = x$purity,
    ref = x$reference_genome
  )
  
  return(x)
}