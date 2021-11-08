#
#



#' CNAqc-based purity-optimisation pipeline for Sequenza CNA calling.
#'
#' @description This pipeline implements optimised allele-specific
#' copy number calling with Sequenza and CNAqc, optimising purity estimated
#' from data, and allele-specific segments.
#'
#' Requirements:
#'
#' 1. (Prerequisites) Run Sequenza utils steps outside R
#'    a. Process a FASTA file to produce a GC Wiggle track file
#'    b. Process BAM and Wiggle files to produce a seqz file:
#'    c. Post-process by binning the original seqz file:
#'
#'    Instructions to customize these steps are available at
#'    \url{https://cran.r-project.org/web/packages/sequenza/vignettes/sequenza.html}
#'
#' 2. Then this pipeline (\code{Sequenza_CNAqc}) will optimize CNA calling and purity
#' estiamtion using peak detection from CNAqc. A maximum number of attempts
#' \code{max_runs} is attempted, starting from a broad set of initial conditions
#' for purity (5% to 100%). Further refinements steps will restrict the Sequenza purity
#' range based on the score computed by \code{analyse_peaks}, which determines a purity
#' gradient. Each run is saved in a folder named \code{run-1}, \code{run-2}, etc.
#' If at any step a \code{PASS} is obtained via \code{analyse_peaks},
#' the pipeline stops and a symbolic link (\code{final}) pointing to the correct
#' Sequenza solution is created. If no run reaches a \code{PASS}  status, the
#' \code{final} folder is never created but the best reported solution is print
#' out to screen. Every CNAqc output object containing somatic mutations called by
#' Sequenza and allele-specific CNAs is available in the specific running folders.
#'
#' @note  Ploidy is kept constant by this pipeline, using a possible range of ploidy values
#' that are specified by the user. If at any time a purity proposal made by CNAqc leads
#' to inconsistent values (i.e., outside $[0,1]$) the pipeline stops prompting to the user
#' to adjust manually the range of \code{ploidy}.
#'
#' @param sample_id The id of the sample.
#' @param seqz_file
#' @param sex
#' @param cellularity
#' @param ploidy
#' @param max_runs
#' @param reference
#' @param normalization.method
#' @param window
#' @param gamma
#' @param kmin
#' @param min.reads.baf
#' @param min.reads
#' @param min.reads.normal
#' @param max.mut.types
#' @param ... Optional parameters passed to the \code{analyse_peaks} function
#' by CNAqc. Tune these to change error tolerance or karyotypes to use for QC.
#'
#' @return Nothing. If the pipeline could fit a correct purity value then
#' a \code{final} folder links to the final Sequenza results. Otherwise the best
#' solution - not final - is reported to screen. If the solution in general can
#' not be found error logs and screen outputs suggest the user how to check,
#' manually, how to improve the fits.
#'
#' @export
#'
#' @examples
#' \donotrun{
#'
#' # Make some example based on https://cran.r-project.org/web/packages/sequenza/vignettes/sequenza.html
#' install.packages("sequenza")
#'
#' # Outside R
#' # pip install sequenza-utils
#'
#' # Process a FASTA file to produce a GC Wiggle track file
#' # sequenza−utils gc_wiggle −w 50 --fasta hg19.fa -o hg19.gc50Base.wig.gz
#'
#' # Process BAM and Wiggle files to produce a seqz file
#' # sequenza−utils bam2seqz -n normal.bam -t tumor.bam --fasta hg19.fa \
#' # -gc hg19.gc50Base.wig.gz -o out.seqz.gz
#'
#' # Post-process by binning the original seqz file:
#' sequenza−utils seqz_binning --seqz out.seqz.gz -w 50 -o out small.seqz.gz
#'
#' # Run this pipeline
#' Sequenza_CNAqc(
#'    sample_id = 'tumour',
#'    seqz_file = 'small.seqz.gz', # Binned file
#'    sex = "F" # if female
#' )
#' }
Sequenza_CNAqc = function(
  sample_id,
  seqz_file,
  sex,
  cellularity = c(0.05, 1),
  ploidy = c(1.8, 5.4),
  max_runs = 3,
  reference = 'GRCh38',
  normalization.method = 'median',
  window = 1e5,
  gamma = 280,
  kmin = 300,
  min.reads.baf = 50,
  min.reads = 50,
  min.reads.normal = 15,
  max.mut.types = 1,
  ...
)
{
  # Auxiliary function: Sequenza outputs parser
  parse_Sequenza_CNAs = function(out, x)
  {
    # Files required
    segments_file = paste0(out, '/', x, "_segments.txt")
    mutations_file = paste0(out, '/', x, "_mutations.txt")
    purity_file = paste0(out, '/', x, "_confints_CP.txt")
    alternatives_file = paste0(out, '/', x, "_alternative_solutions.txt")

    if(!file.exists(segments_file))
      stop(paste0("Missing required segments file: ", segments_file))

    if(!file.exists(mutations_file))
      stop(paste0("Missing required mutations file: ", mutations_file))

    if(!file.exists(purity_file))
      stop(paste0("Missing required solutions file: ", purity_file))

    if(!file.exists(alternatives_file))
      stop(paste0("Missing alternative solutions file: ", alternatives_file))

    # Segments data -- with CNAq format conversion which is used in evoverse
    segments = readr::read_tsv(segments_file, col_types = readr::cols()) %>%
      dplyr::rename(
        chr = chromosome,
        from = start.pos,
        to = end.pos,
        Major = A,
        minor = B
      ) %>%
      dplyr::select(chr, from, to, Major, minor, dplyr::everything())

    # Mutations data -- with CNAq format conversion which is used in evoverse
    mutations = readr::read_tsv(mutations_file, col_types = readr::cols()) %>%
      dplyr::rename(
        chr = chromosome,
        from = position,
      ) %>%
      dplyr::rowwise() %>%
      tidyr::separate(mutation, sep = '>', into = c('ref', 'alt'), remove = FALSE) %>%
      dplyr::mutate(to = from + nchar(alt)) %>%
      dplyr::select(chr, from, to, ref, alt, dplyr::everything())

    # Solution data
    solutions = readr::read_tsv(purity_file, col_types = readr::cols())

    purity = solutions$cellularity[2]
    ploidy = solutions$ploidy.estimate[2]

    # Alternative solutions data
    alternative_solutions = readr::read_tsv(alternatives_file, col_types = readr::cols())

    fits = list()
    fits$input = x
    fits$out = out
    fits$segments = segments
    fits$mutations = mutations
    fits$purity = purity
    fits$ploidy = ploidy
    fits$alternative_solutions = alternative_solutions

    return(fits)
  }

  # Auxiliary function: Sequenza check input parameters
  Sequenza_check_inputs = function(sample_id, seqzFile, sex, cellularity,  ploidy, max_runs, reference)
  {
    if(!is.character(sample_id))
      cli::cli_abort("Unrecogniseable sample id, will not proceed.")

    cli::cli_alert("Sample id: {.field {sample_id}}")

    if(!file.exists(seqzFile))
      cli::cli_abort("File {.field {seqzFile}} does not exist, will not proceed.")

    sz = format(object.size(seqzFile), units = "auto")
    cli::cli_alert("Sequenza seqz file: {.value \"{seqzFile}\"} [{.field {sz}}]")

    if(!(sex %in% c("M", "F")))
      cli::cli_abort("Sex must be either female (F) or male (M).")

    sx = dplyr::case_when(
      sex == "F" ~ "Female (F)",
      sex == "M" ~ "Male (M)",
    )
    cli::cli_alert("Sample sex: {.field {sx}}")

    if(
      (cellularity %>% length() != 2) |
      !(cellularity[1] %>% is.numeric()) |
      !(cellularity[2] %>% is.numeric())
    )
      cli::cli_abort("Cellularity needs to be a numeric vector with two values.")

    cli::cli_alert("Testing cellularity (%) in [{.field {cellularity[1]}}, {.field {cellularity[2]}}]")

    if(
      (ploidy %>% length() != 2) |
      !(ploidy[1] %>% is.numeric()) |
      !(ploidy[2] %>% is.numeric())
    )
      cli::cli_abort("Ploidy needs to be a numeric vector with two values.")

    cli::cli_alert("Testing ploidy (real-valued) in [{.field {ploidy[1]}}, {.field {ploidy[2]}}]")

    if(!(max_runs %>% is.numeric()))
      cli::cli_abort("Maximum number of runs should be numeric.")

    cli::cli_alert("Testing maximum {.field {max_runs}} purity adjustments with CNAqc")

    if(
      !(reference %>% is.character()) |
      !(reference %in% c("GRCh38", "GRCh37", "hg19"))
    )
      cli::cli_abort("Reference must be any of GRCh38 or hg19/GRCh37.")

    cli::cli_alert("Reference genome: {.field {reference}}")
  }

  #### Check for Sequenza library ####
  if(!require(sequenza)) {
    cli::cli_abort("The R Sequenza package is not installed, will not proceed.")
  }

  cli::cli_h1("Sequenza CNA calling wrapper with CNAqc")

  Sequenza_check_inputs(
    sample_id = sample_id,
    seqzFile = seqz_file,
    sex = sex,
    cellularity = cellularity,
    ploidy = ploidy,
    max_runs = max_runs,
    reference = reference
    )

  #### If male use chrY, if female don't ####
  
  is_female <-  sex == "F"
  chromosomes <- paste0("chr", c(1:22, "X"))

  if (!as.logical(is_female)) {
    chromosomes <- c(chromosomes,'chrY')
  }

  #### Parameters that will be logged as RDS ####
  run_params = list(
    sample_id = sample_id,
    seqz_file = seqz_file,
    sex = sex,
    chromosome.list = chromosomes,
    normalization.method = normalization.method,
    window = window,
    gamma = gamma,
    kmin = kmin,
    min.reads.baf = min.reads.baf,
    min.reads = min.reads,
    min.reads.normal = min.reads.normal,
    max.mut.types = max.mut.types,
    max_runs = max_runs
  )

  dir.create('R_logs')
  saveRDS(run_params, 'R_logs/logs.rds')

  # cli::cli_h2("Extraction of required information [{crayon::blue('sequenza.extract')}]")
  
  cat("\n")
  cli::cli_process_start("Extraction of required information [{crayon::blue('sequenza.extract')}]")

  # Run sequenza.extract - with parameters optimised for Genomics England data (100x tumour WGS + 30X normal WGS)
  
  seqzExt <- sequenza::sequenza.extract(
    file = seqz_file,
    chromosome.list = chromosomes,
    normalization.method = normalization.method,
    window = window,
    gamma = gamma,
    kmin = kmin,
    min.reads.baf = min.reads.baf,
    min.reads = min.reads,
    min.reads.normal = min.reads.normal,
    max.mut.types = max.mut.types)

  cli::cli_process_done()

  # Select only the more reliable autosomes for fitting cellularity and ploidy parameters
  chr.fit <- paste0("chr",c(1:22))

  # Perform max_runs attempts
  QC_PASS = FALSE
  out_dir = NULL
  scores = NULL
  cnaqc_objs = NULL

  for(i in 1:max_runs)
  {
    #### Fit cellularity and ploidy parameters - this is where CNAqc is involved ####
    out_dir = paste0('run-', i)
    dir.create(out_dir)

    cat("\n")
    cli::cli_process_start("Run {.field {i}} of fit [{crayon::blue('sequenza.fit')}]")

    # Run sequenza.fit
    paraSpace <- sequenza::sequenza.fit(
      sequenza.extract = seqzExt,
      cellularity = seq(cellularity[1], cellularity[2], 0.01),
      ploidy = seq(ploidy[1], ploidy[2], 0.1),
      chromosome.list = chr.fit,
      female = as.logical(is_female)
      )

    # Dump results
    sequenza::sequenza.results(
      sequenza.extract = seqzExt,
      cp.table = paraSpace,
      sample.id = sample_id,
      out.dir = out_dir,
      female = as.logical(is_female)
      )

    cli::cli_process_done()

    # Parse outputs
    cat("\n")
    cli::cli_process_start("Quality control with CNAqc")

    fits = parse_Sequenza_CNAs(out_dir, x = sample_id)
    
    fits$mutations = fits$mutations %>% rename(VAF="F", DP="good.reads") %>% 
      mutate(NV=as.integer(VAF*DP))

    # Perform QC, and save RDS
    cnaqc_obj = CNAqc::init(
      snvs = fits$mutations,
      cna = fits$segments,
      purity = fits$purity,
      ref = reference
    ) %>%
      CNAqc::analyze_peaks()

    scores = append(scores, cnaqc_obj$peaks_analysis$score)
    cnaqc_objs = append(cnaqc_objs, cnaqc_obj)

    # saveRDS(
    #   file = paste0(out_dir, '/', sample_id, '_cnaqc.rds')
    # )

    # QC status
    col_msg = ifelse(cnaqc_obj$peaks_analysis$QC == "PASS", 'green', 'brown')
    txt_msg = ifelse(cnaqc_obj$peaks_analysis$QC == "PASS",
    "CNAqc QC PASS, will stop running Sequenza",
    "CNAqc QC FAIL, needs to re-run Sequenza")

    cat("\n")
    cat(cli::boxx(
      txt_msg,
      padding = 1,
      float = 'center',
      background_col = col_msg
    ))
    cat("\n")

    print(cnaqc_obj)
    cli::cli_process_done()

    if(cnaqc_obj$peaks_analysis$QC == "PASS")
    {
      QC_PASS = TRUE
      break
    }
    else
    {
      # Purity adjustments
      proposed_purity = fits$purity
      adjusted_purity = (proposed_purity + cnaqc_obj$peaks_analysis$score) %>% round(2)

      # Special case, adjusted purity is incompatible with reality
      # so we stop, raise an error, and ask to change ploidy
      if(adjusted_purity < 0 | adjusted_purity > 1)
      {
        cli::cli_alert_danger("Proposed purity {.value {adjusted_purity}} is outside \\
                       range [0,1] and cannot be used. Please consider \\
                       adjusting input ploidy and re-running.")
        break
      }

      old_cellularity = cellularity %>% round(2)
      cellularity = c(
        adjusted_purity - 0.05/i,
        adjusted_purity + 0.05/i
      ) %>% round(2)

      cli::cli_alert(
        "Cellularity was [{.value {old_cellularity[1]}}, {.value {old_cellularity[2]}}] \\
        adjusted to [{.value {cellularity[1]}}, {.value {cellularity[2]}}] \\
        -> CNAqc prediction {.field {adjusted_purity*100}%}."
      )
    }
  }

  # Summary
  if(QC_PASS)
  {
    cli::cli_alert_success(
      "Joint Seqeuenza/CNAqc analysis succesfull. The symbolic link  \\
      {crayon::blue('final')} will point to the best result obtained so."
    )

    cnaqc_objs[[length(cnaqc_objs)]] %>% print

    R.utils::createLink(link = 'final', target = out_dir)
  }
  else
  {
    cli::cli_alert_danger(
      "Joint Seqeuenza/CNAqc analysis unsuccesfull. Consider re-running
      manually this sample. \n \\
      The best result obtained so far is."
    )

    best_fit = which.min(scores)
    cnaqc_objs[[best_fit]] %>% print
  }
}


