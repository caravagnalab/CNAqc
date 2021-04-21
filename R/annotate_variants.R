#' Annotate variant in a CNAqc object
#'
#' The function provides an easy and fast way to annotate SNVs, in particular we annotate the locations of the different mutations, the consequences for
#' those in coding regions and we compare them to a set of known drivers.
#'
#' @param x a data.frame with coloumns c("chr", "from", "to").
#' @param driver_list a one coloumn data.frame with gene symbols
#' @param polyphen whete to annotate polyphen scores
#' @param filter_tumor_type filter intogen drivers by tumor type
#'
#' @return a data.frame with additiona coloumns `location` with the position of the variant in the genome (coding, intron, threeUTR, ...),
#' `consequence` with the consequence of coding mutations (synonymous, nonsynonymous, ...), `is_driver` a bool that indicate if the gene is a driver,
#' `gene_symbol` for the annotated corresponding gene symbol (if the variant is in a gene) and `driver_label` with the driver label written as
#' gene_refAA->varAA (is NA in case `is_driver` is FALSE).
#'
#' The `chr` coloumn int the input should contain the chr prefix before the chromosome number. The annotation process is based on the package
#' \href{https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html}{VariantAnnotation} check it for more information.
#'
#' @export
#'
#' @examples
#'
#'\dontrun{library(CNAqc)
#'
#'data('example_dataset_CNAqc', package = 'CNAqc')
#'
#' snvs <- example_dataset_CNAqc$snvs
#'
#'snvs_annotated <- annotate_variants(snvs)}
#'

annotate_variants <- function(x,
                              ref = "hg19",
                              driver_list = CNAqc::intogen_drivers,
                              polyphen = FALSE,
                              filter_tumor_type = NULL) {
  tx_pkg <-  paste0("TxDb.Hsapiens.UCSC.", ref, ".knownGene")
  bs_pkg <- paste0("BSgenome.Hsapiens.UCSC.", ref)


  if (!require(tx_pkg,  character.only = T, quietly = T) %>% suppressWarnings()) {
    stop(
      paste0(
        "Please install ",
        tx_pkg,
        " from Bioconducator to use the variant annotation function."
      )
    )
  }

  if (!require(bs_pkg,  character.only = T, quietly = T) %>% suppressWarnings()) {
    stop(
      paste0(
        "Please install ",
        bs_pkg,
        " from Bioconducator to use the variant annotation function."
      )
    )
  }

  if (!require("Organism.dplyr",
               character.only = T,
               quietly = T) %>% suppressWarnings()) {
    stop(
      paste0(
        "Please install ",
        "Organism.dplyr",
        " from Bioconducator to use the variant annotation function."
      )
    )
  }

  if (!require("org.Hs.eg.db",
               character.only = T,
               quietly = T) %>% suppressWarnings()) {
    stop(
      paste0(
        "Please install ",
        "org.Hs.eg.db",
        " from Bioconducator to use the variant annotation function."
      )
    )
  }

  inp = NULL

  if (inherits(x, 'cnaqc'))
  {
    cli::cli_alert("Preparing annotations for the CNAqc object")
    inp = x$snvs
  }
  else
  {
    cli::cli_alert("Preparing annotations for a dataframe")
    inp = x
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
  cli::cli_process_start("Locaring variants with {crayon::yellow('VariantAnnotation')}")

  loc <-
    VariantAnnotation::locateVariants(rd, txdb, VariantAnnotation::AllVariants())

  cli::cli_process_done()

  # Src organism DB
  cli::cli_process_start("Sourcing organism")

  src <- src_organism(tx_pkg)

  cli::cli_process_done()

  cli::cli_process_start("Traslating Entrez ids")

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
    dplyr::select(chr, from, to, consequence, refAA, varAA) %>% unique() %>% group_by(chr, from, to) %>%
    summarize(
      consequence = paste(consequence, collapse = ":"),
      refAA = paste(refAA, collapse = ":"),
      varAA = paste(varAA, collapse = ":")
    ) %>%  ungroup()

  res <-
    dplyr::left_join(loc_df, output_coding,  by = c("chr", "from", "to")) %>% mutate(
      driver_label = paste0(gene_symbol, "_", refAA, "->", varAA)
    )
  res <-
    dplyr::left_join(res,
                     driver_list %>% dplyr::select(gene_symbol, is_driver),
                     by = "gene_symbol")

  res <-  res %>%
    mutate(is_driver = ifelse(
      is.na(is_driver) |
        !grepl(location, pattern = "coding") |
        is.na(consequence) |
        grepl(consequence, pattern = "^(?!non)synonymous", perl = T),
      FALSE,
      TRUE
    )) %>%
    mutate(driver_label = ifelse(is_driver, driver_label, NA)) %>%  unique()

  cli::cli_h3("Found {.green {sum(res$is_driver)}} driver(s)")
  print(res %>% filter(is_driver))


  if (polyphen) {
    # ahaha
    stop("Sorry not yet implemented, wait just a couple of days (I swear)!")
    poly_pkg <- "PolyPhen.Hsapiens.dbSNP131"
    if (!require(poly_pkg,  character.only = T, quietly = T)) {
      stop(
        paste0(
          "Please install ",
          poly_pkg,
          " from Bioconducator to use the PolyPhen variant annotation function."
        )
      )
    }
  }

  final_table = left_join(inp %>% mutate(to = to + 1),
                          res %>% mutate(to = to + 1), by = c("chr", "from", "to")) %>% dplyr::arrange(-is_driver)


  final_table %>% filter(is_driver)
  res %>% filter(is_driver) %>% mutate(to = to + 1)
  inp %>% filter(Hugo_Symbol =='VHL')%>% mutate(to = to + 1) %>% dplyr::select(chr, from, to, ref, alt)

  if (inherits(x, 'cnaqc'))
  {
    x$snvs = final_table
    final_table = x
  }


  return(final_table)
}
