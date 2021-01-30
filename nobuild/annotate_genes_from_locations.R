# This function uses
# GenomicRanges
# GenomicFeatures
# TxDb.Hsapiens.UCSC.hg19.knownGene
# org.Hs.eg.db
annotate_genes_from_locations = function(x)
{
  if (!requireNamespace("Homo.sapiens", quietly = TRUE)) {
    stop(
      "Package \"Homo.sapiens\" needed for this function to work. Please install it.",
      call. = FALSE
    )
  }
  
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop(
      "Package \"GenomicRanges\" needed for this function to work. Please install it.",
      call. = FALSE
    )
  }
  
  if (!requireNamespace("GenomicFeatures", quietly = TRUE)) {
    stop(
      "Package \"GenomicFeatures\" needed for this function to work. Please install it.",
      call. = FALSE
    )
  }
  
  if (!requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE)) {
    stop(
      "Package \"TxDb.Hsapiens.UCSC.hg19.knownGene\" needed for this function to work. Please install it.",
      call. = FALSE
    )
  }
  
  if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
    stop(
      "Package \"TxDb.Hsapiens.UCSC.hg38.knownGene\" needed for this function to work. Please install it.",
      call. = FALSE
    )
  }
  
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop(
      "Package \"org.Hs.eg.db\" needed for this function to work. Please install it.",
      call. = FALSE
    )
  }
  
  if (!all(c('chr', 'from', 'to', 'ref', 'alt') %in% colnames(x$data)))
    stop(
      "Missing genomic coordinates (chr, from, ref, alt) from the input data, cannot annotate genes."
    )
  
  mutations = x$snvs %>%
    dplyr::select(chr, from, ref, alt) %>%
    dplyr::rename(chrom = chr, start = from) %>%
    dplyr::mutate(start = as.numeric(start),
                  end = start + nchar(alt)) %>%
    dplyr::select(chrom, start, end)
  
  GR_df_mutations = GenomicRanges::makeGRangesFromDataFrame(mutations %>%
                                                              data.frame(stringsAsFactors = FALSE))
  
  if (x$reference_genome %in% c("hg38", "GRCh38"))
    genes_list = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  
  if (x$reference_genome %in% c("hg19", "GRCh37"))
    genes_list = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  
  genes_map = org.Hs.eg.db::org.Hs.egSYMBOL
  
  GR_df_overlaps = IRanges::subsetByOverlaps(GenomicFeatures::genes(genes_list),
                                             GR_df_mutations) %>%
    data.frame(stringsAsFactors = FALSE) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(seqnames = paste(seqnames),
                  strand = paste(strand)) %>%
    dplyr::left_join(as.data.frame(genes_map), by = 'gene_id') %>%
    dplyr::select(-strand, -gene_id, -width)
  
  colnames(GR_df_overlaps) = c('chr', 'from', 'to', 'gene')
  colnames(mutations) = c('chr', 'from', 'to')
  
  final_mapping = lapply(1:nrow(GR_df_overlaps),
                         function(x) {
                           mutations %>%
                             dplyr::filter(chr == GR_df_overlaps$chr[x],
                                           from >= GR_df_overlaps$from[x],
                                           to <= GR_df_overlaps$to[x]) %>%
                             dplyr::mutate(gene = GR_df_overlaps$gene[x])
                         }) %>%
    Reduce(f = bind_rows)
  
  mutations_enriched = mutations %>%
    dplyr::left_join(final_mapping, by = c('chr', 'from', 'to')) %>%
    dplyr::distinct(chr, from, to, .keep_all = T)
  
  # table(mutations_enriched$gene) %>% sort(decreasing = TRUE) %>% pioDisp()
  mobster::Clusters(x)  %>%
    dplyr::mutate(from = as.numeric(from),
                  to = from + nchar(ref)) %>%
    dplyr::left_join(mutations_enriched, by = c('chr', 'from', 'to')) %>%
    dplyr::select(chr, from, to, ref, alt, gene, everything())
}