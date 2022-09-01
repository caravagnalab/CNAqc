#' Annotate variant with VariantAnnotation.
#'
#' @description  The function provides an easy and fast way to annotate variants,
#' reporting the locations of the different mutations and the consequences for
#' those in coding regions. Then putative drivers are returned upon matching
#' against lists of known drivers.
#'
#' @param x A data.frame with columns c("chr", "from", "to"), or a CNAqc object.
#' @param driver_list a one column data.frame with gene symbols.
#' @param collapse whether or not to merge different location, consequences and genes in a single string or split the mutations
#' @param filter_tumor_type filter Intogen drivers by tumor type.
#'
#' @return a data.frame with an additional columns `location` reporting the position of the variant in the genome (`coding`, `intron`, `threeUTR`, ...),
#' `consequence` with the consequence of coding mutations (`synonymous`, `nonsynonymous`, ...), `is_driver` a boolean that indicates if the gene is a driver,
#' `gene_symbol` for the annotated corresponding gene symbol (if the variant is in a gene) and `driver_label` with the driver label written as
#' `gene_refAA->varAA` (is `NA` in case `is_driver = FALSE`).
#'
#' The `chr` column in the input should contain the chr prefix before the chromosome number. The annotation process is based on the package
#' \href{https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html}{VariantAnnotation}.
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
                              ref = "hg19",
                              driver_list = CNAqc::intogen_drivers,
                              collapse = TRUE,
                              filter_tumor_type = NULL)
{
  tx_pkg <-  paste0("TxDb.Hsapiens.UCSC.", ref, ".knownGene")
  bs_pkg <- paste0("BSgenome.Hsapiens.UCSC.", ref)

  # Check for available packages
  rqp = function(x) {
    if (!require(x,  character.only = T, quietly = T) %>% suppressWarnings()) {
      cli::cli_abort("Bioconducator package {.field {x}} is required to annotate variants")
    }
  }

  # Required packages
  tx_pkg %>% rqp()
  bs_pkg %>% rqp()
  "Organism.dplyr" %>% rqp()
  "org.Hs.eg.db" %>% rqp()

  # Define input
  inp = x

  if (inherits(x, 'cnaqc'))
  {
    cli::cli_alert_info("Extracting mutations from a CNAqc object.")
    inp = x %>% Mutations()
  }

  # Cancel other annotations
  if(
    any(c('is_driver', "driver_label") %in% colnames(inp))
  )
  {
    cli::cli_alert_warning("There are existing annotations in your data, they will be cancelled.")
    print(inp %>% filter(is_driver))

    inp$is_driver = inp$driver_label = NULL
  }

  # Get data and ranges
  txdb <- eval(parse(text = tx_pkg))
  inp <- inp %>% mutate(to = to - 1)

  rd <-
    GenomicRanges::makeGRangesFromDataFrame(
      inp,
      start.field = "from",
      end.field = "to",
      seqnames.field = "chr",
      keep.extra.columns = F
    )

  # VariantAnnotation package
  cli::cli_process_start("Locating variants with {crayon::yellow('VariantAnnotation')}")

  loc <-
    VariantAnnotation::locateVariants(rd, txdb, VariantAnnotation::AllVariants())

  cli::cli_process_done()

  # Src organism DB
  cli::cli_process_start("Traslating Entrez ids")

  src <- src_organism(tx_pkg)

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
    as_tibble() %>%
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
    group_by(chr, from, to, gene_symbol) %>%
    summarize(location = paste(location, collapse = ":"),
              .groups = 'keep') %>%
    ungroup()

  input_coding <-
    dplyr::left_join(loc_df %>%  filter(grepl(location, pattern = "coding")),
                     inp,
                     by = c("chr", "from", "to"))

  # inp %>% filter(chr == 'chr1', from == 985357)

  input_coding <-  dplyr::left_join(
    input_coding,
    seqlengths(Hsapiens) %>%  as.data.frame() %>%
      tibble::rownames_to_column() %>% dplyr::rename(chr = "rowname", length = "."),
    by = "chr"
  ) %>%
    filter(from < length)

  input_coding_grange <-
    GenomicRanges::makeGRangesFromDataFrame(
      input_coding,
      start.field = "from",
      end.field = "to",
      seqnames.field = "chr",
      keep.extra.columns = F
    )

  cli::cli_process_done()

  cli::cli_process_start("Predicting coding")

  coding <-
    VariantAnnotation::predictCoding(input_coding_grange,
                                     txdb,
                                     seqSource = Hsapiens,
                                     DNAStringSetList(lapply(input_coding$alt, function(i)
                                       i)) %>% unlist)
  cli::cli_process_done()

  if (!is.null(filter_tumor_type)) {
    driver_list <-
      driver_list %>% filter(CANCER_TYPE %in% filter_tumor_type)
  }

  colnames(driver_list)[1] <-  "gene_symbol"
  driver_list$is_driver <-  TRUE
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
    dplyr::select(chr, from, to, consequence, refAA, varAA) %>% unique()

  res <-
    dplyr::left_join(loc_df, output_coding,  by = c("chr", "from", "to")) %>% mutate(
      driver_label = paste0(gene_symbol, "_", refAA, "->", varAA)
    )
  res <-
    dplyr::left_join(res,
                     driver_list %>% dplyr::select(gene_symbol, is_driver),
                     by = "gene_symbol")

  res <-  res %>%
    mutate(is_driver =
             ifelse(
               is.na(is_driver) |
                 !grepl(location, pattern = "coding") |
                 is.na(consequence) |
                 grepl(consequence, pattern = "^(?!non)synonymous", perl = T),
      FALSE,
      TRUE)
    ) %>%
    mutate(driver_label = ifelse(is_driver, driver_label, NA)) %>%  unique()

  cli::cli_h3("Found {.green {sum(res$is_driver)}} driver(s)")
  print(res %>% filter(is_driver))


  if(collapse){
    res <- res %>%  group_by(chr,from,to) %>%
      summarize(gene_symbol = paste(unique(gene_symbol), collapse = ":"),
                location = paste(unique(location), collapse = ":"),
                refAA = paste(unique(refAA), collapse = ":"),
                varAA = paste(unique(varAA), collapse = ":"),
                consequence = paste(unique(consequence), collapse = ":"),
                driver_label= paste(unique(driver_label), collapse = ":"),
                is_driver = any(is_driver))
    res$driver_label <- gsub(res$driver_label, pattern = "NA:", replacement = "", fixed = TRUE)
  }


  final_table = left_join(
    inp %>% mutate(to = to + 1),
    res %>% mutate(to = to + 1),
    by = c("chr", "from", "to")
    ) %>%
    dplyr::arrange(-is_driver)



  # final_table %>% filter(is_driver)
  # res %>% filter(is_driver) %>% mutate(to = to + 1)
  # inp %>% filter(Hugo_Symbol =='VHL')%>% mutate(to = to + 1) %>% dplyr::select(chr, from, to, ref, alt)

  if (inherits(x, 'cnaqc'))
  {
    x$mutations = final_table
    final_table = x
  }

  return(final_table)
}
