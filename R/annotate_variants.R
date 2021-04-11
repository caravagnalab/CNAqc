

annotate_variant <- function(x, driver_list = CNAqc::intogen_drivers){

  ref <-  x$reference
  tx_pkg <-  paste0("TxDb.Hsapiens.UCSC.",ref, ".knownGene")



  if(!require(tx_pkg,  character.only = T, quietly = T)){
    stop(paste0("Please install ", tx_pkg, " from Bioconducator to use the variant annotation function."))
  }

  if(!require(tx_pkg,  character.only = T, quietly = T)){
    stop(paste0("Please install ", tx_pkg, " from Bioconducator to use the variant annotation function."))
  }

  if(!require("Organism.dplyr",  character.only = T, quietly = T)){
    stop(paste0("Please install ", "Organism.dplyr", " from Bioconducator to use the variant annotation function."))
  }

  if(!require("org.Hs.eg.db",  character.only = T, quietly = T)){
    stop(paste0("Please install ", "org.Hs.eg.db", " from Bioconducator to use the variant annotation function."))
  }

  library(tx_pkg, character.only = T)

  txdb <- eval(parse(text = tx_pkg))

  rd <-  GenomicRanges::makeGRangesFromDataFrame(x$snvs, start.field = "from", end.field = "to", seqnames.field = "chr", keep.extra.columns = F, ignore.strand = T)

  loc <- VariantAnnotation::locateVariants(rd, txdb, VariantAnnotation::AllVariants(),
                                           ignore.strand = T)

  loc_df <- as.data.frame(loc) %>%
    dplyr::mutate(segment_id = paste(seqnames, start, end, sep = ":"), chr = seqnames, from = start, to = end, location = LOCATION) %>%
    dplyr::select(chr, from, to, location, GENEID)

  head(loc_df)

  src <- src_organism(tx_pkg)

  AnnotationDbi::select(src, keys=loc_df$GENEID, columns=c("symbol", "tx_name"), keytype="entrez")

}



