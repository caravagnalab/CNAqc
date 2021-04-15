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
  # print(to_find)

  found = drivers_list %>%
    filter(driver_gene %in% to_find) %>%
    pull(driver_gene) %>%
    unique

  cli::cli_alert_info("{length(found)} found in the input list of {length(drivers_list$driver_gene) %>% unique()} genes.")
  # if(length(found) > 0) print(found)

  if(length(found) == 0) return(x)

  which_driver = intersect(
    which(x$snvs[[gene_column]] %in% drivers_list$driver_gene),
    which(x$snvs[[function_column]] %in% function_entries)
  )

  x$snvs$is_driver = FALSE
  x$snvs$is_driver[which_driver] = TRUE

  cli::cli_h3('Drivers')
  x$snvs %>%
    filter(is_driver) %>%
    print

  # cli::cli_alert_info("Reassembling a new CNAqc object with these drivers")
  # x = CNAqc::init(
  #   x$snvs,
  #   cna = x$cna,
  #   purity = x$purity,
  #   ref = x$reference_genome
  # )

  # cli::cli_alert_info("Updating the CNAqc object with these drivers")
  return(x)
}


# exones_GRCh38

# exones = exones_hg19
#
# exones_idxs = easypar::run(
#   FUN = function(i)
#   {
#     w = exones %>%
#       dplyr::filter(chr == x$snvs$chr[i],
#                     from <= x$snvs$from[i],
#                     to >= x$snvs$to[i]) %>%
#       distinct(gene)
#
#     if (nrow(w) > 0) {
#       # cli::cli_alert_success("Found exonic {.field {w$gene %>% paste(collapse = ':')}} mutation.")
#       return(data.frame(idx = i , ensemble_genes = paste(w$gene, collapse = ':')))
#     }
#     else
#       return(NULL)
#   },
#   PARAMS = lapply(1:nrow(x$snvs), list),
#   parallel = FALSE,
#   filter_errors = F,
#   export = 'exones'
# )
# exones_idxs = Reduce(bind_rows, exones_idxs)
#
# x$snvs$ensembl_exon = FALSE
# x$snvs$ensembl_exon[exones_idxs$idx] = TRUE
# x$snvs$ensembl_exon_gene = NA
# x$snvs$ensembl_exon_gene[exones_idxs$idx] = exones_idxs$ensemble_genes
#
# x$snvs %>% filter(ensembl_exon) %>% dplyr::select(GENE, ANNOVAR_FUNCTION, starts_with("ensembl_exon"))
#
# drv_list = CNAqc::DriverDBv3 %>%
#   distinct(driver_gene)
#
#   rename(ensembl_exon_gene = driver_gene)
#
# x$snvs = x$snvs %>%
#   left_join(drv_list, by = 'ensembl_exon_gene') %>%
#   mutate(is_driver = !is.na(n_tools))
#
# x$snvs %>% filter(is_driver) %>% dplyr::select(GENE, ANNOVAR_FUNCTION, starts_with("ensembl_exon"), cancer, n_tools) %>% View



# exones_GRCh38 %>% filter(gene =='TP53') %>% pull(from) %>% min
# exones_GRCh38 %>% filter(gene =='TP53') %>% pull(from) %>% max
# gene_coordinates_GRCh38 %>% filter(gene =='TP53')




