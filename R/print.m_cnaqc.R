#' Print for class \code{'m_cnaqc'}.
#'
#' @param x An obj of class \code{'m_cnaqc'}.
#' @param ... Default S3 method parameter.
#'
#' @return Nothing.
#'
#' @exportS3Method print m_cnaqc
#' @export print.m_cnaqc
#'
#' @examples
#' 
#' ### search for an example to add here and include in the package
#' 
#' \dontrun{
#' data('example_multisample', package = 'CNAqc')
#'
#' print(example_multisample)
#' }
#' 

print.m_cnaqc <- function(x, ...) {
  
  stopifnot(inherits(x, "m_cnaqc"))
  
  n_segments_clonal = lapply(x, function(y) y$mutations_on_shared$cna %>% nrow) %>% unlist %>% unique()
  n_segments_subclonal = lapply(x, function(y) y$mutations_on_shared$cna_subclonal %>% nrow) %>% unlist %>% unique()
  tot_segments = n_segments_subclonal + n_segments_clonal
  ref_gen = lapply(x, function(y) y$mutations_on_shared$reference_genome) %>% unlist %>% unique()
  
  
  n_mutations_shared = lapply(x, function(y) y$mutations_on_shared$mutations %>% nrow) %>% unlist
  n_mutations_private = lapply(x, function(y) y$mutations_on_private %>% nrow) %>% unlist

  cli::cli_rule(paste(crayon::bgCyan(crayon::white("[ multi CNAqc ]")), 
                      "{.field {tot_segments}} total shared segments ({.field {n_segments_clonal}} clonal, {.field {n_segments_subclonal}}) subclonal. Rerefence genome: {.field {ref_gen}}"))
  
  # cli::cli_ul(paste(names(old_segments), old_segments, sep = " = "))
  
  cli::cli_h3("Clonal mutations")
  cli::cli_ul(paste(names(n_mutations_private), paste(crayon::green(n_mutations_shared), "mutations on shared segments and", crayon::green(n_mutations_private), "mutations on private segments"), sep = ": "))
  
  cli::cli_h3("Avaiable information and analyses")
  cli::cli_ul(lapply(x, function(y) {
    y$original_additional_info %>% names
  }) %>% unlist() %>% unique)
}
