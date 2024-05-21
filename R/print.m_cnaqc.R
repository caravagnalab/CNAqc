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
#########################################################################################

# S3 method print for class obj m_cnaqc

print.m_cnaqc <- function(x, ...) {
  stopifnot(inherits(x, "m_cnaqc"))
  
  n_segments_total = x$m_cnaqc_stats$n_cna_new_segmentation %>% unique
  n_mutations_total = x$m_cnaqc_stats$n_mutations_new_segmentation %>% unique
  
  # n_segments_clonal = lapply(x, function(y) y$mutations_on_shared$cna %>% nrow) %>% unlist %>% unique()
  # n_segments_subclonal = lapply(x, function(y) y$mutations_on_shared$cna_subclonal %>% nrow) %>% unlist %>% unique()
  # tot_segments = n_segments_subclonal + n_segments_clonal
  ref_gen = lapply(x$cnaqc_obj_new_segmentation, function(y)
    y$reference_genome) %>% unlist %>% unique()
  samples = get_sample_name(x) %>% paste(., collapse = ", ")
  
  # n_mutations_shared = lapply(x, function(y) y$mutations_on_shared$mutations %>% nrow) %>% unlist
  # n_mutations_private = lapply(x, function(y) y$mutations_on_private %>% nrow) %>% unlist
  
  cli::cli_rule(crayon::bgCyan(crayon::white("[ multi CNAqc ]")))
  cat("\n")
  cli::cli_ul(
    c(
      "Rerefence genome: {.field {ref_gen}}",
      "Samples: {.field {samples} }",
      "{.field {n_segments_total}} total shared segments across samples",
      "{.field {n_mutations_total}} mutations on shared positions across samples"
      
    )
  )
  
  
  # cli::cli_ul(paste(names(old_segments), old_segments, sep = " = "))
  
  # cli::cli_h3("Clonal mutations")
  # cli::cli_ul(paste(names(n_mutations_shared), paste(crayon::green(n_mutations_shared), "mutations on shared positions"), sep = ": "))
  
  # lapply(names(x), function(y) {
  if ("original_cnaqc_objc" %in% names(x)) {
    cli::cli_h3("Original CNAqc object avaiable for {.field {samples}}")
  } else {
    cli::cli_h3("Original CNAqc object not avaiable for {.field {samples}}")
  }
  # %>% invisible()
}

######################################################
#' Plot for class \code{'m_cnaqc'}.
#'
#' @description
#'
#' The default plot is the CNA segments in wide format
#'
#' @param x An obj of class \code{'m_cnaqc'}.
#' @param ... Default S3 method parameter.
#'
#' @return Nothing.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' 
#' }
#' 
############################################################################
# plot S3 method for m_cnaqc obj

plot.m_cnaqc <- function(x, ...)
{
  stopifnot(inherits(x, "m_cnaqc"))
  
  shared = get_sample(x, sample = get_sample_name(x), which_obj = "shared")
  CNAqc::plot_multisample_CNA(shared, layout = "circular")
}  
