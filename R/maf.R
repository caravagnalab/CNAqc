#' Import MAF annotations.
#'
#' @description
#'
#' This function imports
#' \href{https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/#introduction}{Mutation Annotation Formats (MAF)}
#' data into a CNAqc object, or to a mutations
#' dataframe in the format ready for CNAqc. Therefore the expected uses are
#'
#' \itemize{
#'  \item{"[Option 1]"}{Mutations + CNA + Purity -> CNAqc -> MAF annotation}
#'  \item{"[Option 2]"}{Mutations -> MAF annotation -> Mutations (with MAF) + CNA + Purity -> CNAqc}
#' }
#'
#' MAF annotations should have been created using the
#' \href{https://github.com/mskcc/vcf2maf}{vcf2maf}
#' utility. At this point, if one has created the input data for CNAqc from the
#' original VCF file, the MAF can be added to the CNAqc object (Option 1), or to
#' the available mutations before creating the CNAqc obejct (Option 2).
#'
#' Mutations are associated based on genome locations and substitutions. Since
#' MAFs tend to use some different convention for genome location (especially
#' for indels), mutations are matched by closes genomic coordinates.
#'
#' Note that if one has many CNAqc objects with their MAF annotations, a cohort of
#' MAFs can be exported and functions from Bioconductor package
#' \href{https://bioconductor.org/packages/release/bioc/html/maftools.html}{maftools}
#' can be used to plot data from multiple patients.
#'
#' @param x A CNAqc object, or a dataframe of mutations in the input format
#' for CNAqc.
#' @param maf The file MAF associated containing the annotations in MAF format.
#'
#' @seealso function \code{\link{as_maftools_cohort}} to convert multiple CNAqc
#' objects with MAF annotations into a single MAF cohort;  package
#' \href{https://bioconductor.org/packages/release/bioc/html/maftools.html}{maftools}
#' to summarize, analyze and visualize MAF Files; the utility
#' \href{https://github.com/mskcc/vcf2maf}{vcf2maf} to create MAF files from VCFs,
#' using the
#' \href{https://www.ensembl.org/info/docs/tools/vep/index.html}{Ensembl Variant Effect Predictor (VEP)}
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
augment_with_maf = function(x, maf)
{
  # Manage inputs
  if(inherits(x, 'cnaqc'))
  {
    cli::cli_h1("Augmenting a CNAqc object with its MAF")

    if(x %>% has_MAF_annotations())
    {
      cli::cli_alert_warning("Input has already MAF annotations that will be dropped")

      cn = colnames(x$mutations)
      cn = grepl("MAF.", cn)

      x$mutations = x$mutations[, !cn, drop = FALSE]
    }

    cli::cli_h2("Input CNAqc object")

    x %>% print()
  }

  if(x %>% is.data.frame())
  {
    cli::cli_h1("Augmenting a mutation dataframe with MAF annotations")
    cat("\n")

    required_columns = c('chr', 'from', 'to', 'ref', 'alt', 'DP', 'NV', 'VAF')

    if (!all(required_columns %in% colnames(x)))
      stop("Bad mutation format, see the manual.")

    print(x)
  }

  # Mutations table
  mutations_table = NULL

  if(inherits(x, 'cnaqc'))
    mutations_table = x %>% Mutations() # so we take those for subclonal CNAs as well

  if(x %>% is.data.frame())
    mutations_table = x

  # MAF loading
  cli::cli_h2("MAF input")

  MAF_input = maftools::read.maf(maf)

  cat("\n")
  print(MAF_input)

  # MAF conversion
  MAF_input = MAF_input@data %>% tibble::as_tibble()

  shared_colnames = intersect(mutations_table %>% colnames,
                              MAF_input %>% colnames)

  cn = colnames(MAF_input)
  cn_match = !(cn %in% shared_colnames)
  colnames(MAF_input)[cn_match] = paste0("MAF.", cn[cn_match])

  # MAF_input = MAF_input %>%
  #   dplyr::mutate(
  #     chr = MAF.Chromosome,
  #     from = MAF.Start_Position,
  #     to = MAF.End_Position + 1,
  #     # re-scaling required
  #     ref = MAF.Reference_Allele,
  #     alt = ifelse(
  #       MAF.Tumor_Seq_Allele1 == ref,
  #       MAF.Tumor_Seq_Allele2,
  #       MAF.Tumor_Seq_Allele1
  #     )
  #   )

  # Matching in MAF style
  MAF_input$maf_match_id = MAF_input$match_distance_maf = NA
  mutations_table$maf_match_id = 1:nrow(mutations_table)

  for(i in 1:nrow(MAF_input))
  {
    f_cnaqc = mutations_table %>%
      dplyr::filter(chr == MAF_input$MAF.Chromosome[i]) %>%
      dplyr::mutate(dist_maf = abs(from - MAF_input$MAF.Start_Position[i])) %>%
      dplyr::arrange(dist_maf)

    # MAF variant i-th matches CNAqc variant f_cnaqc$id[1]
    # with distance f_cnaqc$dist_maf[1]
    MAF_input$maf_match_id[i] = f_cnaqc$maf_match_id[1]
    MAF_input$match_distance_maf[i] = f_cnaqc$dist_maf[1]
  }

  if(any(is.na(MAF_input$maf_match_id)))
  {
    cli::cli_alert_warning("Some MAF mutations are not matched")
    MAF_input %>%
      dplyr::filter(is.na(maf_match_id)) %>%
      dplyr::select(
        MAF.Hugo_Symbol,
        MAF.Chromosome,
        MAF.Start_Position,
        MAF.End_Position,
        MAF.Variant_Type,
        MAF.Variant_Classification
      ) %>%
      print()
  }

  if(any(MAF_input$match_distance_maf > 0)) {
    cli::cli_alert_warning("The following MAF mutations are assigned, but not matched exactly -- might be indels?")

    MAF_input %>%
      dplyr::filter(match_distance_maf > 0) %>%
      dplyr::arrange(dplyr::desc(match_distance_maf)) %>%
      dplyr::select(
        MAF.Hugo_Symbol,
        MAF.Chromosome,
        MAF.Start_Position,
        MAF.End_Position,
        MAF.Variant_Type,
        MAF.Variant_Classification
        ) %>%
      print()
  }

  mutations_table =
    mutations_table %>%
    dplyr::left_join(MAF_input, by = c("maf_match_id", shared_colnames))

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

  return(x)
}



#' Convert a CNAqc object to a maftools object.
#'
#' @description
#'
#' A single CNAqc object for which MAF annotations have been added using function
#' \code{\link{augment_with_maf}}, can be converted as an object for the Bioconductor
#' package [maftools](https://bioconductor.org/packages/release/bioc/html/maftools.html).
#'
#' The export will use only driver mutations data or all the annotated mutations,
#' depending on the input flag `only_drivers`. Moreover, CNAqc copy number states
#' of a list of desired genes can also be augmented to the
#' [maftools](https://bioconductor.org/packages/release/bioc/html/maftools.html)
#' object.
#'
#' @param x A CNAqc object with MAF annotations.
#' @param only_drivers If `TRUE`, only driver mutations are used, otherwised all.
#' When `TRUE`, if drivers are not annotated, an error is thrown.
#' @param CNA_genes The list of genes for ....
#' @param assembly If `TRUE`...
#'
#' @return
#'
#' @seealso function \code{\link{augment_with_maf}} to add MAF annotations to a
#' CNAqc object, to be used before running `as_maftools_obj`.
#'
#' @export
#'
#' @examples
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
#'    # Extraction
#'    x %>% as_maftools_obj
#' }
as_maftools_obj = function(x,
                           only_drivers = TRUE,
                           CNA_genes = NULL,
                           cross_reference = TRUE,
                           assembly = TRUE)
{
  if (!(inherits(x, 'cnaqc'))) {
    cli::cli_abort("The input object is not a CNAqc object.")
  }

  if (!(x %>% has_MAF_annotations)) {
    cli::cli_abort("MAF annotations are missing, cannot use this function.")
  }

  cli::cli_h1("Conversion to maftools")

  mutations_data = x %>% Mutations()
  colnames(mutations_data) = gsub('MAF.', '', mutations_data %>% colnames)

  if (only_drivers)
  {
    if (!(x %>% has_drivers)) {
      cli::cli_abort("Driver data is missing, cannot use only_drivers = TRUE.")
    }

    mutations_data = mutations_data %>% dplyr::filter(is_driver)
    cli::cli_alert("Using {.field {nrow(mutations_data)}} driver mutations")

  } else
    cli::cli_alert("Using {.field {nrow(mutations_data)}} mutations")

  if (!is.null(CNA_genes))
  {
    cli::cli_alert("Extracting CNA data for {.field {length(CNA_genes)}} genes")

    if(cross_reference)
      CNA_genes = intersect(CNA_genes, mutations_data$Hugo_Symbol)

    # CNA mapping
    genes_map = CNA_gene(x, genes = CNA_genes)

    genes_map = genes_map %>%
      dplyr::rename(Gene = gene,
                    CN = karyotype) %>%
      dplyr::mutate(Sample_name = x$sample) %>%
      dplyr::select(-chr, -from, -to, -Major, -minor) %>%
      dplyr::select(Gene, Sample_name, CN)

    if (assembly)
    {
      return(maftools::read.maf(mutations_data,
                                cnTable = genes_map,
                                verbose = FALSE))
    }
    else
    {
      return(list(mutations = mutations_data,
                  CNA = genes_map))
    }

  }

  # No CNA data
  if (assembly)
  {
    return(maftools::read.maf(mutations_data,
                              verbose = FALSE))
  }
  else
  {
    return(list(mutations = mutations_data, CNA = NULL))
  }
}


#' Convert a list of CNAqc object to a maftools object.
#'
#' @description
#'
#' Like \code{\link{as_maftools_obj}} but for multiple CNAqc objects, it creates
#' a unique [maftools](https://bioconductor.org/packages/release/bioc/html/maftools.html)
#' object with all the samples at once.
#'
#' Parameters have the same meaning as for function \code{\link{as_maftools_obj}}.
#'
#' @param x A list of CNAqc objects with MAF annotations.
#' @param only_drivers If `TRUE`, only driver mutations are used, otherwised all.
#' When `TRUE`, if drivers are not annotated, an error is thrown.
#' @param CNA_genes
#'
#' @return
#'
#' @seealso function \code{\link{augment_with_maf}} to add MAF annotations to a
#' CNAqc object, to be used before running `as_maftools_obj`.
#'
#' @export
#'
#' @examples
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
#'    # Extraction
#'    x %>% as_maftools_obj
#' }
as_maftools_cohort = function(x,
                              only_drivers = TRUE,
                              CNA_genes = NULL)
{
  inputs = lapply(x,
                  as_maftools_obj,
                  only_drivers = only_drivers,
                  CNA_genes = CNA_genes,
                  cross_reference = TRUE,
                  assembly = FALSE)

  pooled_mutations = lapply(inputs, function(x) x[[1]]) %>% Reduce(f = dplyr::bind_rows)
  pooled_CNA = lapply(inputs, function(x) x[[2]]) %>% Reduce(f = dplyr::bind_rows)

  pooled_CNA = pooled_CNA %>% mutate(CN = ifelse(CN %in% c("1:0", "2:0", "2:1", "2:2"), CN, 'Other'))


  if(!is.null(pooled_CNA) & nrow(pooled_CNA) > 0)
  {
    return(maftools::read.maf(pooled_mutations,
                              cnTable = pooled_CNA,
                              verbose = FALSE))
  }

  return(maftools::read.maf(pooled_mutations,
                            verbose = FALSE))

  # lapply(x,
  #        as_maftools_obj,
  #        only_drivers = only_drivers,
  #        CNA_genes = CNA_genes) %>%
  #   maftools::merge_mafs()
}

#' Extract per-gene copy number status.
#'
#' @description
#'
#' By default, this function maps a list of genes to their copy number status
#' using clonal CNA segments. The genes used are all the known human genes,
#' whose coordinates are available inside the CNAqc package and are identified
#' by common names (e.g., TP53). The function can restrict to a subset of genes
#' (faster computation) if one passes, via the `genes` parameter, a vector of
#' gene symbols.
#'
#' @param x A CNAqc object.
#' @param genes Optional, a vector of gene symbols of interest. If `NULL`, all
#' the human genes are used, according to the genome reference of input `x`.
#'
#' @return A tibble with columns `gene` (gene name), `from`/`to` (gene delimiters),
#' `Major`/`minor`/`karyotype` as the information for the copy number segment
#' where the gene sits. Note that if the gene maps to a subclonal segment
#' this is not returned.
#'
#' @export
#'
#' @examples
#' # Example input data released with the package
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' print(example_dataset_CNAqc)
#'
#' # Note the outputs to screen
#' x = init(mutations = example_dataset_CNAqc$mutations, cna = example_dataset_CNAqc$cna, purity = example_dataset_CNAqc$purity)
#'
#' # Get mapping for all the known human genes - takes a bit longer
#' CNA_gene(x)
#'
#' # Use known genes set
#' CNA_gene(x, genes = c("APC", "KRAS", "NRAS", "TP53"))
CNA_gene = function(x, genes = NULL)
{
  # Coordinates
  if (x$reference_genome %in% c("hg38", "GRCh38"))
  {
    data('gene_coordinates_GRCh38', package = 'CNAqc')
    coordinates = gene_coordinates_GRCh38
  }
  else{
    if (x$reference_genome %in% c("hg19", "GRCh37"))
    {
      data('gene_coordinates_hg19', package = 'CNAqc')
      coordinates = gene_coordinates_hg19
    }
    else
      cli::cli_abort("Unrecognised reference -- cannot extract gene copy numbers")
  }

  if (!is.null(genes))
    coordinates = coordinates %>% dplyr::filter(gene %in% genes)

  # Clonal CNAs
  cna = x %>% CNA(type = 'clonal')

  # Scan CNA by gene or viceversa depends on what is faster
  if (nrow(cna) < length(genes))
  {
    # cli::cli_alert("Scanning by  ")

    genes_map = easypar::run(
      FUN = function(i) {
        coordinates %>%
          dplyr::filter(chr == cna$chr[i], from >= cna$from[i], to <= cna$to[i]) %>%
          dplyr::mutate(
            Major = cna$Major[i],
            minor = cna$minor[i],
            karyotype = paste0(Major, ":", minor)
          )
      },
      PARAMS = lapply(1:nrow(cna), list),
      parallel = FALSE
    )

    genes_map = Reduce(dplyr::bind_rows, genes_map)

    return(genes_map)
  }
  else
  {
    genes_map = easypar::run(
      FUN = function(i) {
        where = cna %>%
          dplyr::filter(chr == coordinates$chr[i],
                        from <= coordinates$from[i],
                        to >= coordinates$to[i])

        data.frame(
          i = i,
          Major = where$Major[1],
          minor = where$minor[1],
          karyotype =
            paste0(where$Major[1], ":", where$minor[1])
        )
      },
      PARAMS = lapply(1:length(genes), list),
      parallel = FALSE
    )

    genes_map = Reduce(dplyr::bind_rows, genes_map)

    # what found
    coordinates = coordinates %>%
      dplyr::mutate(i = dplyr::row_number()) %>%
      dplyr::left_join(genes_map, by = 'i') %>%
      dplyr::select(-i)

    return(coordinates)
  }
}
