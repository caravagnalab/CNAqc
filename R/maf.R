#' Import MAF annotations into a CNAqc object
#'
#' @description
#'
#' This function import MAF annotations into a CNAqc object. MAF annotations
#' should have been created using the [vcf2maf](https://github.com/mskcc/vcf2maf)
#' utility. At this point, if one has created the input data for CNAqc from the
#' original VCF file, the MAF can be added to the CNAqc object. Mutations are
#' associated based on genome locations and substitutions. Once many CNAqc objects
#' have been augmented with their MAF annotations, a cohort of MAFs can be
#' exported and functions from Bioconductor package [maftools](https://bioconductor.org/packages/release/bioc/html/maftools.html)
#' can be used to plot data from multiple patients.
#'
#' @param x A CNAqc object.
#' @param maf The file MAF associated containing the annotations in MAF format
#' for the input object `x`.
#'
#' @seealso function \code{\link{as_maftools_cohort}} to convert multiple CNAqc
#' objects with MAF annotations into a single MAF cohort;  package
#' [maftools](https://bioconductor.org/packages/release/bioc/html/maftools.html) to
#' summarize, analyze and visualize MAF Files; the utility [vcf2maf](https://github.com/mskcc/vcf2maf)
#' to create MAF files from VCFs, using the [Ensembl Variant Effect Predictor (VEP)](https://www.ensembl.org/info/docs/tools/vep/index.html)
#' utility.
#'
#' @return A CNAqc object like `x` where the mutations are associated to the
#' MAF annotations, matched by genomic coordinates of the mutations. The S3
#' print method for `x` will report the presence of MAF annotaitons.
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
#'    # check they are in (there should be many columns with "MAF." prefix)
#'    x %>% Mutations %>% colnames
#' }
augment_with_maf = function(x, maf)
{
  cli::cli_h1("Augmenting a CNAqc object with its MAF")

  cli::cli_h2("Input CNAqc object")

  x %>% print()

  # MAF loading
  cli::cli_h2("MAF input")

  MAF_input = maftools::read.maf(maf)

  cat("\n")
  print(MAF_input)

  # MAF conversion
  MAF_input = MAF_input@data %>% as_tibble()

  shared_colnames = intersect(x$mutations %>% colnames,
                              MAF_input %>% colnames)

  cn = colnames(MAF_input)
  colnames(MAF_input)[!(cn %in% shared_colnames)] = paste0("MAF.", cn[!(cn %in% shared_colnames)])

  MAF_input = MAF_input %>%
    dplyr::mutate(
      chr = MAF.Chromosome,
      from = MAF.Start_Position,
      to = MAF.End_Position + 1,
      # re-scaling required
      ref = MAF.Reference_Allele,
      alt = ifelse(
        MAF.Tumor_Seq_Allele1 == ref,
        MAF.Tumor_Seq_Allele2,
        MAF.Tumor_Seq_Allele1
      )
    )


  x$mutations = x$mutations %>%
    dplyr::left_join(MAF_input,
                     by = c("chr", "from", "to", "ref", "alt", shared_colnames))

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
as_maftools_obj = function(x,
                           only_drivers = TRUE,
                           CNA_genes = NULL)
{

  if(!(inherits(x, 'cnaqc'))){
    cli::cli_abort("The input object is not a CNAqc object.")
  }

  if(!(x %>% has_MAF_annotations)){
    cli::cli_abort("MAF annotations are missing, cannot use this function.")
  }

  cli::cli_h1("Conversion to maftools")

  mutations_data = x %>% Mutations()
  colnames(mutations_data) = gsub('MAF.', '', mutations_data %>% colnames)

  if (only_drivers)
  {
    if(!(x %>% has_driver_data)){
      cli::cli_abort("Driver data is missing, cannot use only_drivers = TRUE.")
    }

    mutations_data = mutations_data %>% dplyr::filter(is_driver)
    cli::cli_alert("Using {.field {nrow(mutations_data)}} driver mutations")

  } else
    cli::cli_alert("Using {.field {nrow(mutations_data)}} mutations")

  if (!is.null(CNA_genes))
  {
    cli::cli_alert("Extracting CNA data for {.field {length(CNA_genes)}} genes")

    # CNA mapping
    genes_map = CNA_gene(x, genes = CNA_genes)

    genes_map = genes_map %>%
      dplyr::rename(Gene = gene,
                    CN = karyotype) %>%
      dplyr::mutate(Sample_name = x$sample) %>%
      dplyr::select(-chr,-from,-to,-Major,-minor) %>%
      dplyr::select(Gene, Sample_name, CN)

    m =  maftools::read.maf(mutations_data,
                            cnTable = genes_map,
                            verbose = FALSE)
  }
  else
    m =  maftools::read.maf(mutations_data,
                            verbose = FALSE)

  return(m)
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
  lapply(x,
         as_maftools_obj,
         only_drivers = only_drivers,
         CNA_genes = CNA_genes) %>%
    maftools::merge_mafs()
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
  if(nrow(cna) < length(genes))
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
