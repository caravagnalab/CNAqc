#' Annotate variants and drivers.
#'
#' @description  Easy and fast way to annotate input mutations, and detect
#' potential driver mutations. This function computes the locations of
#' the different mutations and the consequences of substituions mapped to
#' coding regions, using \href{https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html}{VariantAnnotation} and other Bioconductor packages.
#' Then, putative drivers are annotated upon matching from an input list that,
#' by default, is compiled from the \href{https://www.intogen.org}{Intogen} database. Drivers are selected among
#' coding substituions with known effect.
#'
#' @param x A CNAqc object.
#' @param drivers A dataframe in the format of the `intogen_drivers` one released with
#' @param make_0_span Remove -1 from the coloumn `to`, in case SNPs is indicated by `from` and `to = from + 1` to make
#' it effectively a point mutation.
#' @param collapse if the same mutation has more than one consqeunce or location in different transcript, it collapse them
#' by concatenating them using `:`
#' CNAqc. In particular, it must contain a column named `gene` to identify gene names.
#'
#' @return A CNAqc object with variants annotated. For each variant this object contains:
#'
#' - `location` reporting the position of the variant in the genome (`coding`, `intron`, `threeUTR`, ...),
#' - `consequence` with the consequence of coding mutations (`synonymous`, `nonsynonymous`, ...),
#' - `is_driver` a boolean that indicates if the gene is a driver,
#' - `gene_symbol` for the annotated corresponding gene symbol (if the variant is in a gene)
#' - `driver_label` with the driver label written as `gene_refAA->varAA` (`NA` in case `is_driver = FALSE`).
#'
#' The annotation process is based on the package \href{https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html}{VariantAnnotation}.
#'
#' @export
#'
#' @examples
#'
#'\dontrun{
#' library(CNAqc)
#'
#' data('example_dataset_CNAqc', package = 'CNAqc')
#'
#' mutations <- example_dataset_CNAqc$mutations
#'
#' mutations_annotated <- annotate_variants(mutations)
#' }
#'
annotate_variants <- function(x,
                              drivers = CNAqc::intogen_drivers,
                              make_0_span = FALSE,
                              collapse = TRUE)
{
  # Old params, now fixed.
  collapse = TRUE

  inputs = annotate_variants_preprocess(x)
  mutations = inputs$mutations
  reference = inputs$reference

  if(reference == "GRCh38") reference = "hg38"
  if(reference == "GRCh37") reference = "hg19"

  # Check for available packages
  rqp = function(x) {
    if (!require(x,  character.only = T, quietly = T) %>% suppressWarnings()) {
      cli::cli_abort("Bioconducator package {.field {x}} is required to annotate variants")
    }
  }

  # Required packages
  tx_pkg <-  paste0("TxDb.Hsapiens.UCSC.", reference, ".knownGene")
  bs_pkg <- paste0("BSgenome.Hsapiens.UCSC.", reference)

  tx_pkg %>% rqp()
  bs_pkg %>% rqp()
  "Organism.dplyr" %>% rqp()
  "org.Hs.eg.db" %>% rqp()
  "GenomicRanges" %>% rqp()

  # Get data and ranges
  txdb <- eval(parse(text = tx_pkg))

  if(make_0_span) mutations <- mutations %>% mutate(from = to - 1)

  rd <-
    GenomicRanges::makeGRangesFromDataFrame(
      mutations,
      start.field = "from",
      end.field = "to",
      seqnames.field = "chr",
      keep.extra.columns = F
    )
  rd <- rd[as.character(seqnames(rd)) %in% unique(as.character(seqlevels(txdb)))]
  seqlevels(rd) <- as.character(unique(seqnames(rd)))
  
  
  # VariantAnnotation package
  cli::cli_process_start("Locating variants with {crayon::yellow('VariantAnnotation')}")

  loc <-
    VariantAnnotation::locateVariants(rd, txdb, VariantAnnotation::AllVariants())

  cli::cli_process_done()

  # Src organism DB
  cli::cli_process_start("Traslating Entrez ids")

  src <- src_organism(tx_pkg, overwrite = TRUE)

  translated_enttrz <-
    AnnotationDbi::select(
      src,
      keys = loc$GENEID,
      columns = c("symbol"),
      keytype = "entrez"
    )

  cli::cli_process_done()

  cli::cli_process_start("Transforming data")

  map <-
    setNames(translated_enttrz$entrez, object =  translated_enttrz$symbol)

  loc$gene_symbol <- map[loc$GENEID]

  loc_df <- as.data.frame(loc) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      segment_id = paste(seqnames, start, end, sep = ":"),
      chr = seqnames,
      from = start,
      to = end,
      location = LOCATION
    ) %>%
    dplyr::select(chr, from, to, location, gene_symbol) %>%
    dplyr::filter(gene_symbol != "NA") %>%
    unique() %>%
    dplyr::group_by(chr, from, to, gene_symbol) %>%
    dplyr::summarize(location = paste(location, collapse = ":"),
              .groups = 'keep') %>%
    dplyr::ungroup()

  mutationsut_coding <-
    dplyr::left_join(loc_df %>%
                       dplyr::filter(grepl(location, pattern = "coding")),
                     mutations %>% mutate(from = as.integer(from), to = as.integer(to)),
                     by = c("chr", "from", "to"))

  mutationsut_coding <-  dplyr::left_join(
    mutationsut_coding,
    seqlengths(Hsapiens) %>%
      as.data.frame() %>%
      tibble::rownames_to_column() %>%
      dplyr::rename(chr = "rowname", length = "."),
    by = "chr"
  ) %>%
    filter(from < length)

  mutationsut_coding_grange <-
    GenomicRanges::makeGRangesFromDataFrame(
      mutationsut_coding,
      start.field = "from",
      end.field = "to",
      seqnames.field = "chr",
      keep.extra.columns = F
    )

  cli::cli_process_done()

  cli::cli_process_start("Predicting coding")

  coding <-
    VariantAnnotation::predictCoding(mutationsut_coding_grange,
                                     txdb,
                                     seqSource = Hsapiens,
                                     DNAStringSetList(lapply(mutationsut_coding$alt, function(i)
                                       i)) %>% unlist)

  output_coding <- as.data.frame(coding) %>%
    dplyr::mutate(
      segment_id = paste(seqnames, start, end, sep = ":"),
      chr = seqnames,
      from = start,
      to = end,
      consequence = CONSEQUENCE,
      refAA = REFAA,
      varAA = VARAA
    ) %>%
    dplyr::select(chr, from, to, consequence, refAA, varAA) %>%
    unique()

  cli::cli_h3("Coding substitutions found")
  print(output_coding)
  cat("\n")

  cli::cli_process_done()

  # Drivers annotation
  cli::cli_process_start("Drivers annotation")

  if(!is.data.frame(drivers)) stop("`drivers` has to be a dataframe ")

  driver_genes = drivers$gene %>% unique

  res <-
    dplyr::left_join(loc_df, output_coding,  by = c("chr", "from", "to")) %>%
    dplyr::mutate(
      driver_label = paste0(gene_symbol, "_", refAA, "->", varAA),
      is_driver = gene_symbol %in% driver_genes
    ) %>%
    dplyr:: mutate(is_driver =
             ifelse(
               is.na(is_driver) |
                 !grepl(location, pattern = "coding") |
                 is.na(consequence) |
                 grepl(consequence, pattern = "^(?!non)synonymous", perl = T),
               FALSE,
               gene_symbol %in% driver_genes # check against the list
               )
    ) %>%
    mutate(driver_label = ifelse(is_driver, driver_label, NA)) %>%
    unique()

  cli::cli_h3("Found {.green {sum(res$is_driver)}} driver(s)")

  print(res %>% filter(is_driver))

  cli::cli_process_done()

  if(collapse){
    res <- res %>%
      group_by(chr,from,to) %>%
      summarize(gene_symbol = paste(unique(gene_symbol), collapse = ":"),
                location = paste(unique(location), collapse = ":"),
                refAA = paste(unique(refAA), collapse = ":"),
                varAA = paste(unique(varAA), collapse = ":"),
                consequence = paste(unique(consequence), collapse = ":"),
                driver_label= paste(unique(driver_label), collapse = ":"),
                is_driver = any(is_driver)
                ) %>%
      ungroup()

    res$driver_label <- gsub(res$driver_label, pattern = "NA:", replacement = "", fixed = TRUE)
  }

  # Re-assembly CNAqc obect
  final_table = dplyr::left_join(
    mutations,
    res,
    by = c("chr", "from", "to")
  ) %>%
    dplyr::arrange(-is_driver)

  if(make_0_span)  final_table <- final_table %>% dplyr::mutate(to = to + 1)

  x_new = CNAqc::init(
    mutations = final_table,
    cna = x %>% CNA(),
    purity = x$purity,
    ref = x$reference_genome
  )

  return(x_new)
}

annotate_variants_preprocess = function(x)
{
  cli::cli_process_start("Preparing mutations")

  # Define mutationsuts
  mutations = x %>% Mutations()
  reference = x$reference_genome

  # Cancel other annotations
  if(
    any(c('is_driver', "driver_label") %in% colnames(mutations))
  )
  {
    cli::cli_alert_warning("Existing driver annotations in your data will be cancelled.")
    cat("\n")
    print(mutations %>% dplyr::filter(is_driver))
    cat("\n")

    mutations$is_driver = mutations$driver_label = NULL
  }

  cli::cli_process_done()

  return(list(mutations = mutations, reference = reference))
}

# annotate_variants <- function(x,
#                               ref = "hg19",
#                               driver_list = CNAqc::intogen_drivers,
#                               collapse = TRUE,
#                               filter_tumor_type = NULL)
# {
#   tx_pkg <-  paste0("TxDb.Hsapiens.UCSC.", ref, ".knownGene")
#   bs_pkg <- paste0("BSgenome.Hsapiens.UCSC.", ref)
#
#   # Check for available packages
#   rqp = function(x) {
#     if (!require(x,  character.only = T, quietly = T) %>% suppressWarnings()) {
#       cli::cli_abort("Bioconducator package {.field {x}} is required to annotate variants")
#     }
#   }
#
#   # Required packages
#   tx_pkg %>% rqp()
#   bs_pkg %>% rqp()
#   "Organism.dplyr" %>% rqp()
#   "org.Hs.eg.db" %>% rqp()
#
#   # Define mutationsut
#   mutations = x
#
#   if (inherits(x, 'cnaqc'))
#   {
#     cli::cli_alert_info("Extracting mutations from a CNAqc object.")
#     mutations = x %>% Mutations()
#   }
#
#   # Cancel other annotations
#   if(
#     any(c('is_driver', "driver_label") %in% colnames(mutations))
#   )
#   {
#     cli::cli_alert_warning("There are existing annotations in your data, they will be cancelled.")
#     print(mutations %>% filter(is_driver))
#
#     mutations$is_driver = mutations$driver_label = NULL
#   }
#
#   # Get data and ranges
#   txdb <- eval(parse(text = tx_pkg))
#   mutations <- mutations %>% mutate(to = to - 1)
#
#   rd <-
#     GenomicRanges::makeGRangesFromDataFrame(
#       mutations,
#       start.field = "from",
#       end.field = "to",
#       seqnames.field = "chr",
#       keep.extra.columns = F
#     )
#
#   # VariantAnnotation package
#   cli::cli_process_start("Locating variants with {crayon::yellow('VariantAnnotation')}")
#
#   loc <-
#     VariantAnnotation::locateVariants(rd, txdb, VariantAnnotation::AllVariants())
#
#   cli::cli_process_done()
#
#   # Src organism DB
#   cli::cli_process_start("Traslating Entrez ids")
#
#   src <- src_organism(tx_pkg)
#
#   translated_enttrz <-
#     AnnotationDbi::select(
#       src,
#       keys = loc$GENEID,
#       columns = c("symbol"),
#       keytype = "entrez"
#     )
#
#   cli::cli_process_done()
#
#   cli::cli_process_start("Transforming data")
#
#   map <-
#     setNames(translated_enttrz$entrez, object =  translated_enttrz$symbol)
#
#   loc$gene_symbol <- map[loc$GENEID]
#
#   loc_df <- as.data.frame(loc) %>%
#     as_tibble() %>%
#     dplyr::mutate(
#       segment_id = paste(seqnames, start, end, sep = ":"),
#       chr = seqnames,
#       from = start,
#       to = end,
#       location = LOCATION
#     ) %>%
#     dplyr::select(chr, from, to, location, gene_symbol) %>%
#     dplyr::filter(gene_symbol != "NA") %>%
#     unique() %>%
#     group_by(chr, from, to, gene_symbol) %>%
#     summarize(location = paste(location, collapse = ":"),
#               .groups = 'keep') %>%
#     ungroup()
#
#   mutationsut_coding <-
#     dplyr::left_join(loc_df %>%  filter(grepl(location, pattern = "coding")),
#                      mutations,
#                      by = c("chr", "from", "to"))
#
#   # mutations %>% filter(chr == 'chr1', from == 985357)
#
#   mutationsut_coding <-  dplyr::left_join(
#     mutationsut_coding,
#     seqlengths(Hsapiens) %>%  as.data.frame() %>%
#       tibble::rownames_to_column() %>% dplyr::rename(chr = "rowname", length = "."),
#     by = "chr"
#   ) %>%
#     filter(from < length)
#
#   mutationsut_coding_grange <-
#     GenomicRanges::makeGRangesFromDataFrame(
#       mutationsut_coding,
#       start.field = "from",
#       end.field = "to",
#       seqnames.field = "chr",
#       keep.extra.columns = F
#     )
#
#   cli::cli_process_done()
#
#   cli::cli_process_start("Predicting coding")
#
#   coding <-
#     VariantAnnotation::predictCoding(mutationsut_coding_grange,
#                                      txdb,
#                                      seqSource = Hsapiens,
#                                      DNAStringSetList(lapply(mutationsut_coding$alt, function(i)
#                                        i)) %>% unlist)
#   cli::cli_process_done()
#
#   if (!is.null(filter_tumor_type)) {
#     driver_list <-
#       driver_list %>% filter(CANCER_TYPE %in% filter_tumor_type)
#   }
#
#   colnames(driver_list)[1] <-  "gene_symbol"
#   driver_list$is_driver <-  TRUE
#   output_coding <- as.data.frame(coding) %>%
#     dplyr::mutate(
#       segment_id = paste(seqnames, start, end, sep = ":"),
#       chr = seqnames,
#       from = start,
#       to = end,
#       consequence = CONSEQUENCE,
#       refAA = REFAA,
#       varAA = VARAA
#     ) %>%
#     dplyr::select(chr, from, to, consequence, refAA, varAA) %>% unique()
#
#   res <-
#     dplyr::left_join(loc_df, output_coding,  by = c("chr", "from", "to")) %>% mutate(
#       driver_label = paste0(gene_symbol, "_", refAA, "->", varAA)
#     )
#   res <-
#     dplyr::left_join(res,
#                      driver_list %>% dplyr::select(gene_symbol, is_driver),
#                      by = "gene_symbol")
#
#   res <-  res %>%
#     mutate(is_driver =
#              ifelse(
#                is.na(is_driver) |
#                  !grepl(location, pattern = "coding") |
#                  is.na(consequence) |
#                  grepl(consequence, pattern = "^(?!non)synonymous", perl = T),
#       FALSE,
#       TRUE)
#     ) %>%
#     mutate(driver_label = ifelse(is_driver, driver_label, NA)) %>%  unique()
#
#   cli::cli_h3("Found {.green {sum(res$is_driver)}} driver(s)")
#   print(res %>% filter(is_driver))
#
#
#   if(collapse){
#     res <- res %>%  group_by(chr,from,to) %>%
#       summarize(gene_symbol = paste(unique(gene_symbol), collapse = ":"),
#                 location = paste(unique(location), collapse = ":"),
#                 refAA = paste(unique(refAA), collapse = ":"),
#                 varAA = paste(unique(varAA), collapse = ":"),
#                 consequence = paste(unique(consequence), collapse = ":"),
#                 driver_label= paste(unique(driver_label), collapse = ":"),
#                 is_driver = any(is_driver))
#     res$driver_label <- gsub(res$driver_label, pattern = "NA:", replacement = "", fixed = TRUE)
#   }
#
#
#   final_table = left_join(
#     mutations %>% mutate(to = to + 1),
#     res %>% mutate(to = to + 1),
#     by = c("chr", "from", "to")
#     ) %>%
#     dplyr::arrange(-is_driver)
#
#
#
#   # final_table %>% filter(is_driver)
#   # res %>% filter(is_driver) %>% mutate(to = to + 1)
#   # mutations %>% filter(Hugo_Symbol =='VHL')%>% mutate(to = to + 1) %>% dplyr::select(chr, from, to, ref, alt)
#
#   if (inherits(x, 'cnaqc'))
#   {
#     x$mutations = final_table
#     final_table = x
#   }
#
#   return(final_table)
# }

