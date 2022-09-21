#' Import VEP annotations.
#'
#' @description
#'
#' This function imports
#' \href{https://www.ensembl.org/info/docs/tools/vep/index.html}{Ensembl Variant Effect Predictor (VEP)}
#' annotations into a CNAqc object, or to a mutations
#' dataframe in the format ready for CNAqc. Therefore the expected uses are
#'
#' \itemize{
#'  \item{"[Option 1]"}{Mutations + CNA + Purity -> CNAqc -> VEP annotation}
#'  \item{"[Option 2]"}{Mutations -> VEP annotation -> Mutations (with VEP) + CNA + Purity -> CNAqc}
#' }
#'
#' VEP annotations should have been created using the scripts discussed at
#' the
#' \href{https://github.com/caravagnalab/CNAqc}{CNAqc webpage}
#' utility. At this point, if one has created the input data for CNAqc from the
#' original VCF file, the VEP can be added to the CNAqc object (Option 1), or to
#' the available mutations before creating the CNAqc obejct (Option 2).
#'
#' VEP annotations and input mutations are matched by exact genomic coordinates.
#'
#' @param x A CNAqc object, or a dataframe of mutations in the input format
#' for CNAqc.
#' @param maf The file MAF associated containing the annotations in MAF format.
#'
#' @seealso function \code{\link{as_maftools_cohort}} to convert multiple CNAqc
#' objects with MAF annotations into a single MAF cohort;  package
#' [maftools](https://bioconductor.org/packages/release/bioc/html/maftools.html) to
#' summarize, analyze and visualize MAF Files; the utility [vcf2maf](https://github.com/mskcc/vcf2maf)
#' to create MAF files from VCFs, using the [Ensembl Variant Effect Predictor (VEP)](https://www.ensembl.org/info/docs/tools/vep/index.html)
#' utility.
#'
#' @return It depends on `x`
#' \itemize{
#'  \item{"[Option 1]"}{
#'  A CNAqc object like `x` where the mutations are associated to the
#'  MAF annotations, if `x` is a CNAqc object. In this case the S3 print method for `x` will report the
#' presence of MAF annotations.
#' }
#'  \item{"[Option 2]"}{
#'  If `x` is a dataframe, the same dataframe augmented with MAF annotations.
#' }
#'
#' @export
#'
#' @examples
#' # Example with a CNAqc input object
#' if(FALSE)
#' {
#'    # Create your CNAqc object (omissis here) from an original "file.vcf"
#'    x = init(mutations = ..., cna = ..., purity = ...)
#'
#'    # Offline, create your MAF annotations as file "file_vcf.maf" from "file.vcf"
#'    # vcf2maf file.vcf .... file_vcf.maf
#'
#'    # Import into R/CNAqc
#'    x = augment_with_maf(x, maf = "file_vcf.maf")
#'
#'    # check they are in (there should be many columns with "MAF." prefix)
#'    x %>% Mutations %>% colnames
#' }
augment_with_vep = function(x,
                            vep)
{
  # Input selection
  if(inherits(x, 'cnaqc'))
  {
    cli::cli_h1("Augmenting a CNAqc object with VEP annotations")
    cat("\n")

    if(x %>% has_VEP_annotations())
    {
      cli::cli_alert_warning("Input has already VEP annotations that will be dropped")

      cn = colnames(x$mutations)
      cn = grepl("VEP.", cn)

      x$mutations = x$mutations[, !cn, drop = FALSE]
    }

    cli::cli_h2("Input CNAqc object")

    x %>% print()
  }

  if(x %>% is.data.frame())
  {

    cli::cli_h1("Augmenting a mutation dataframe with VEP annotations")
    cat("\n")

    required_columns = c('chr', 'from', 'to', 'ref', 'alt', 'DP', 'NV', 'VAF')

    if (!all(required_columns %in% colnames(x)))
      stop("Bad mutation format, see the manual.")

    print(x)
  }

  # VEP loading
  cli::cli_h2("VEP input")

  VEP_input = readr::read_tsv(
    vep,
    col_names = c(
      "VEP.chr",
      "VEP.from",
      "VEP.ref",
      "VEP.alt",
      "VEP.filter",
      "VEP.annotation",
      "VEP.impact",
      "VEP.gene_name",
      "VEP.gene_code",
      "VEP.exon"
    ),
    col_types = readr::cols()
  ) %>%
    dplyr::mutate(
      VEP.to = VEP.from + nchar(VEP.alt),
      # vep_id_matching = dplyr::row_number()
      ) %>%
    dplyr::select(
      VEP.chr, VEP.from, VEP.to, VEP.ref, VEP.alt, dplyr::everything()
    )

  # cat("\n")
  print(VEP_input)

  # Mutations table
  mutations_table = NULL

  if(inherits(x, 'cnaqc'))
    mutations_table = x %>% Mutations() # so we take those for subclonal CNAs as well

  if(x %>% is.data.frame())
    mutations_table = x

  # VEP can be matched by position
  mutations_table = mutations_table %>%
    dplyr::mutate(
      vep_match_id = paste(chr, from, to, ref, alt, sep = ":")
    ) %>%
    dplyr::left_join(
      VEP_input %>%
        dplyr::mutate(
          vep_match_id = paste(
            VEP.chr,
            VEP.from,
            VEP.to,
            VEP.ref,
            VEP.alt, sep = ":")
        ),
      by = 'vep_match_id'
    )

  # # Closest-matching (as we do for MAFs) -- TOO slow does not make sense
  # VEP_input$vep_match_id = VEP_input$match_distance_vep = NA
  # mutations_table$vep_match_id = 1:nrow(mutations_table)
  #
  # for(i in 1:nrow(VEP_input))
  # {
  #   f_cnaqc = mutations_table %>%
  #     dplyr::filter(chr == VEP_input$VEP.chr[i]) %>%
  #     dplyr::mutate(dist_vep = abs(from - VEP_input$VEP.from[i])) %>%
  #     dplyr::arrange(dist_vep)
  #
  #   # MAF variant i-th matches CNAqc variant f_cnaqc$id[1]
  #   # with distance f_cnaqc$dist_maf[1]
  #   VEP_input$vep_match_id[i] = f_cnaqc$vep_match_id[1]
  #   VEP_input$match_distance_vep[i] = f_cnaqc$dist_vep[1]
  # }

  # Mark exonic mutations
  mutations_table = mutations_table %>%
    dplyr::mutate(
      VEP.EXONIC = ifelse(VEP.exon != "." & !is.na(VEP.exon), TRUE, FALSE)
      )

  # Restore inside the object if required
  if(inherits(x, 'cnaqc'))
  {
    cli::cli_alert_info("Assemblying and returning a new CNAqc object")

    new_objs = init(
      mutations = mutations_table,
      cna = x %>% CNA(),
      purity = x$purity,
      ref = x$reference_genome,
      sample = x$sample
    )

    return(new_objs)
  }

  return(mutations_table)
}
