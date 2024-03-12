#' Augment SBS data for mutational signatures deconvolution
#'
#' @description
#'
#' Augment a CNAqc object with Single Base Substitution (SBS) data for
#' mutational signature deconvolution.
#'
#' This creates set of objects that are stored inside field `SBS` of the input object.
#' These are:
#'
#' * `SNVs`, a tibble for SNVs used for SBS data preparation;
#' * `GRanges`, a GRanges object for the above tibble;
#' * `counts`, a trinucleotide count matrix for SBS deconvolution.
#'
#' In particular, `counts` is the canonnical format to run SBS deconvolution
#' in a variety of tools.
#'
#' @param x A CNAqc object.
#'
#' @return A CNAqc object with required SBS data in the inner field `SBS`.
#' @export
#' @examples
#' \dontrun{
#' data('example_dataset_CNAqc')
#' x = init(mutations = example_dataset_CNAqc$mutations, cna =example_dataset_CNAqc$cna, purity = example_dataset_CNAqc$purity)
#'
#' x = SBS(x)
#'
#' # All SBS data
#' print(x$SBS)
#' }
SBS = function(x)
{
  matrix_counts = function(x){

    data = x

    # create 96 trinucleotides mutation categories
    categories_context <- NULL
    categories_alt <- rep(c(rep("C>A",4),rep("C>G",4),rep("C>T",4),rep("T>A",4),rep("T>C",4),rep("T>G",4)),4)
    categories_cat <- NULL
    cont <- 0
    for(i in c("A","C","G","T")) {
      for(j in 1:6) {
        for(k in c("A","C","G","T")) {
          cont <- cont + 1
          categories_context <- c(categories_context,paste0(k,":",i))
          categories_cat <- c(categories_cat,paste0(k,"[",categories_alt[cont],"]",i))
        }
      }
    }
    mutation_categories <- data.table(context=categories_context,alt=categories_alt,cat=categories_cat)

    # count number of mutations per sample for each category
    data <- merge(mutation_categories[,.(cat)],data.table(sample=data$sample,cat=data$cat)[,.N,by=.(sample,cat)],by="cat",all=TRUE)
    data <- dcast(data,sample~cat,value.var="N")
    data <- data[!is.na(sample),drop=FALSE]
    data[is.na(data)] <- 0

    # make trinucleotides counts matrix
    samples_names <- data$sample
    data <- as.matrix(data[,2:ncol(data),drop=FALSE])
    rownames(data) <- samples_names
    data <- data[sort(rownames(data)),,drop=FALSE]
    data <- data[,sort(colnames(data)),drop=FALSE]
    trinucleotides_counts <- array(0,c(nrow(data),96))
    rownames(trinucleotides_counts) <- rownames(data)
    colnames(trinucleotides_counts) <- sort(as.character(mutation_categories$cat))
    rows_contexts <- rownames(data)
    cols_contexts <- colnames(trinucleotides_counts)[which(colnames(trinucleotides_counts)%in%colnames(data))]
    trinucleotides_counts[rows_contexts,cols_contexts] <- data[rows_contexts,cols_contexts]

    # return trinucleotides counts matrix
    return(trinucleotides_counts)
  }

  stopifnot(inherits(x, "cnaqc"))

  # Load required reference OBJ
  reference = 'none'
  if(x$reference_genome == "GRCh38") reference = "hg38"
  if(x$reference_genome == "GRCh37") reference = "hg19"

  rqp = function(x) {
    if (!require(x,  character.only = T, quietly = T) %>% suppressWarnings()) {
      cli::cli_abort("Bioconductor package {.field {x}} is required to annotate variants")
    }
  }

  bs_pkg <- paste0("BSgenome.Hsapiens.UCSC.", reference)

  cli::cli_alert_info("Using reference package: {.field {bs_pkg}}")

  bs_pkg %>% rqp()
  reference = get(bs_pkg)

  # Extract counts
  sample = x$sample

  # consider only single nucleotide variants involving (A,C,G,T) bases
  snvs = x %>%
    Mutations(type = "SNV") %>%
    dplyr::mutate(sample = x$sample) %>%
    dplyr::select(sample, chr, from, ref, alt) %>%
    dplyr::rename(pos = from) %>%
    as.data.frame()

  # snvs$chr = gsub(snvs$chr, pattern = "chr", replacement = "")
  snvs$id = paste(snvs$chr, snvs$pos, snvs$ref, snvs$alt, sep = ':')

  snvs <- snvs[order(snvs[,"sample"], snvs[,"chr"], snvs[,"pos"]), , drop=FALSE]

  ##### ISSUE 1
  # C'è un motivo per il quale vuoi usare la versione dei reference di UCSC?
  # Tra UCSC e GRCh c'è la differenza della parola chr nei nomi dei cromosomi
  # Per questo esempio ti aggiungo chr nei nomi come fix, però forse andrebbe gestita la cosa
  # Vedi tu come procedere
  # snvs$chr = paste0("chr", snvs$chr)
  #####

  cli::cli_alert_info("Processing {.field n = {nrow(snvs)}} SNVs via GRanges")

  # convert data to GRanges
  data <- GenomicRanges::GRanges(
    as.character(snvs$chr),
    IRanges::IRanges(start = (as.numeric(snvs$pos) - 1), width = 3),
    ref = Biostrings::DNAStringSet(as.character(snvs$ref)),
    alt = Biostrings::DNAStringSet(as.character(snvs$alt)),
    sample = as.character(snvs$sample),
    id = as.character(snvs$id)
  )

  # check that all chromosomes match reference
  if(length(setdiff(seqnames(data), GenomeInfoDb::seqnames(reference))) > 0)
  {
    warning("Check chromosome names, not all match reference genome.")
  }

  # find context for each mutation
  data$context <- getSeq(reference, data)

  # check for any mismatch with BSgenome context
  if(any(subseq(data$context,2,2)!=data$ref)) {
    warning("Check reference bases, not all match context.")
  }

  # get complements and reverse complements
  data$cref <- complement(data$ref)
  data$calt <- complement(data$alt)
  data$rccontext <- reverseComplement(data$context)

  # identify trinucleotides motif
  data$cat <- ifelse(
    data$ref%in%c("C","T"),
    paste0(
      subseq(data$context,1,1), "[", data$ref, ">", data$alt, "]", subseq(data$context,3,3)
    ),
    paste0(
      subseq(data$rccontext,1,1),"[",data$cref,">",data$calt,"]",subseq(data$rccontext,3,3)
    )
  )

  sbs = data %>%
    as_tibble() %>%
    dplyr::select(id,context) %>%
    dplyr::rename(SBS_context = context)

  # Update mutations
  # snvs = snvs %>%
  #   as_tibble() %>%
  #   left_join(sbs, by = 'id')

  ##### ISSUE 2
  ### questa riga da errore perchè non hai campo ID in x$mutations %>% filter(type == "SNV")
  ### ti scrivo fix, poi vedi tu come vuoi sistemarlo
  ###x$sbs = full_join(x$mutations %>% filter(type == "SNV"),sbs, by = "id")
  # tmp = x$mutations %>% filter(type == "SNV")
  # tmp$id = paste0(tmp$chr,":",tmp$from,":",tmp$ref,":",tmp$alt)
  # x$sbs = full_join(tmp,sbs, by = "id")
  # x$GRanges = data
  #####

  snvs = full_join(snvs, sbs, by = "id") %>%
    as_tibble()

  snvs$SBS_substitution = paste0(snvs$ref, '>', snvs$alt)

  schar = function(x) substr(x = x, start  = 1, stop = 1)
  echar = function(x) substr(x = x, start  = nchar(x), stop = nchar(x))

  snvs = snvs %>%
    mutate(
      SBS_label = paste0(schar(SBS_context), '[', SBS_substitution, ']', echar(SBS_context))
    )

  GRanges = data

  x$SBS = list(
    SNVs = snvs,
    GRanges = GRanges,
    counts = matrix_counts(GRanges)
    )

  return(x)
}

#' Extract SBS count data
#'
#' @description
#'
#' Extract SBS data for 96 contexts and 6 possible substitutions, in matrix
#' format.
#'
#' @param x A CNAqc object.
#'
#' @return A one-row matrix with SBS data for 96 contexts and 6 possible substitutions.
#' @export
#' @examples
#' \dontrun{
#' data('example_dataset_CNAqc')
#' x = init(mutations = example_dataset_CNAqc$mutations, cna =example_dataset_CNAqc$cna, purity = example_dataset_CNAqc$purity)
#'
#' x = SBS(x)
#'
#' # SBS analysis
#' print(SBS_counts(x))
#' }
SBS_counts = function(x){

  if(!"SBS" %in% names(x))
  {
    stop("SBS data are not available and should be created - see ?SBS")
  }

  x$SBS$counts
}

#' Plost SBS counts
#'
#' @description Canonical plot of SBS count from a sample, obtained
#' upon mapping substitutions to contexts via function \code{SBS}.
#'
#' The function allows to highlight in red a substitution that occurs with
#' a frequency above a desired cutoff (default 5%).
#'
#' @param x A CNAqc object.
#' @param highlight Minimum frequency to color in red a bar of the barplot (default 5%).
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc')
#' x = init(mutations = example_dataset_CNAqc$mutations, cna =example_dataset_CNAqc$cna, purity = example_dataset_CNAqc$purity)
#'
#' x = SBS(x)
#'
#' # SBS plot
#' plot_SBS(x)
plot_SBS = function(x, highlight = 0.05)
{
  if(!"SBS" %in% names(x)) return(eplot())

  sigs = x$SBS$counts %>%
    reshape2::melt() %>%
    as_tibble() %>%
    dplyr::rename(
      Sample = Var1,
      Feature = Var2,
      Value = value
    )

  # Remove parenthesis
  sigs$substitution = stringr::str_extract_all(sigs$Feature, "\\[[^\\]\\[]*]") %>% unlist()
  sigs$substitution = gsub(pattern = '\\[', replacement = '', x = sigs$substitution)
  sigs$substitution = gsub(pattern = '\\]', replacement = '', x = sigs$substitution)

  # Detect highlights
  sigs = sigs %>%
    dplyr::mutate(
      context = paste0(substr(Feature, 1, 1), '_', substr(Feature, 7, 7)),
      Prop = Value / sum(Value),
      highlight = Prop > highlight
      )

  # Extract a crazy map to colour the reference nucleotide
  sigs$ref = substr(sigs$substitution, 1, 1)

  sigs = sigs %>%
    rowwise() %>%
    mutate(context =
             gsub(
               '_',
               ref,
               context
             ))

    # Nice plot
    lbl = paste0("highlight >", highlight * 100, '%')

    plt =
      ggplot2::ggplot(sigs) +
      ggplot2::geom_hline(
        yintercept = highlight,
        size = .2,
        linetype = 'dashed',
        color = "black"
      )  +
      ggplot2::geom_bar(
        ggplot2::aes(x = context, y = Value, fill = highlight),
        stat="identity",
        position="identity") +
      ggplot2::facet_wrap(~substitution, nrow = 1, scales="free_x") +
        # ggplot2::facet_grid(paste(lbl)~substitution, scales="free", drop = FALSE, space="free_x") +
      my_ggplot_theme() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)
      ) +
      ggplot2::labs(
        y = "SBS count",
        title = paste0(x$sample, ' (n = ', sum(sigs$Value), ')')
      ) +
      ggplot2::scale_fill_manual(values = c(`TRUE` = "indianred3", `FALSE` = 'black')) +
      ggplot2::guides(fill = "none") +
      ggplot2::scale_y_continuous(
        sec.axis = ggplot2::sec_axis(~./sum(sigs$Value),
                                     labels = scales::percent,
                                     name = paste("SBS frequency")
        )
        ) +
      ggplot2::geom_hline(
          yintercept = highlight * sum(sigs$Value),
          color = 'indianred3',
          linetype = 'dashed'
        ) +
      ggplot2::xlab("SBS context")



    return(plt)
}


