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
#' 2. This pipeline (\code{Sequenza_CNAqc}) optimizes CNA calling and purity
#' estimations using peak detection from CNAqc. Starting from a broad set of initial conditions
#' for purity (5\% to 100\%) the pipeline fits cellularity and ploidy values and dumps results.
#' Results are then quality controlled by peak detection via CNAqc. The analysis with CNAqc 
#' can be carried out using either somatic mutations called by Sequenza, or an external set 
#' of input calls. 
#' 
#' 3. At every run the best solution by Sequenza is stored in a cache, and up to two alternative 
#' solutions are generated: one optional by Sequenza, which might propose different exact values for 
#' cellularity and ploidy; one by CNAqc, adjusting the current cellularity of the best Sequenza 
#' solution. The alternative solutions are enqueued in a list L of cellularity and ploidy values 
#' to be tested, based on their presence in the cache. 
#' 
#' 4. Until the list L is empty, the point values of cellularity and ploidy are 
#' tested repeating step 3, using ranges around the proposed values specified by
#' input parameters. Further refinement steps will restrict the Sequenza purity
#' range based on the score computed by \code{analyze_peaks}, which determines a purity
#' gradient. 
#' 
#' 5. Each run is saved in a folder named \code{run-1}, \code{run-2}, etc.
#' When the list L is empty the pipeline stops and a symbolic link \code{final} 
#' pointing to the best Sequenza solution, based on CNAqc score is created. 
#'
#' @note If at any time a purity proposal made by CNAqc leads
#' to inconsistent values (i.e., outside $[0,1]$) the pipeline stops prompting to the user
#' to adjust manually the range of \code{ploidy}.
#'
#' @param sample_id The id of the sample.
#' @param seqz_file A binned \code{.seqz} file.
#' @param mutations If not NULL (default), an external set of mutation calls,
#' in tibble format. Must include colummns "chr", "from", "to", "NV", "NR", "VAF".
#' @param sex Sex of the patient from whom the sample is drawn.
#' @param cellularity Input to \code{sequenza.fit}, see \url{https://cran.r-project.org/web/packages/sequenza/vignettes/sequenza.html}.
#' @param ploidy Input to \code{sequenza.fit}, see \url{https://cran.r-project.org/web/packages/sequenza/vignettes/sequenza.html}. 
#' @param reference Reference genome among supported by CNAqc (GRCh38 or hg19/GRCh37).
#' @param normalization.method Input to \code{sequenza.extract}, see \url{https://cran.r-project.org/web/packages/sequenza/vignettes/sequenza.html}.
#' @param window Input to \code{sequenza.extract}, see \url{https://cran.r-project.org/web/packages/sequenza/vignettes/sequenza.html}.
#' @param gamma Input to \code{sequenza.extract}, see \url{https://cran.r-project.org/web/packages/sequenza/vignettes/sequenza.html}.
#' @param kmin Input to \code{sequenza.extract}, see \url{https://cran.r-project.org/web/packages/sequenza/vignettes/sequenza.html}.
#' @param min.reads.baf Input to \code{sequenza.extract}, see \url{https://cran.r-project.org/web/packages/sequenza/vignettes/sequenza.html}.
#' @param min.reads Input to \code{sequenza.extract}, see \url{https://cran.r-project.org/web/packages/sequenza/vignettes/sequenza.html}.
#' @param min.reads.normal Input to \code{sequenza.extract}, see \url{https://cran.r-project.org/web/packages/sequenza/vignettes/sequenza.html}.
#' @param max.mut.types Input to \code{sequenza.extract}, see \url{https://cran.r-project.org/web/packages/sequenza/vignettes/sequenza.html}.
#' @param delta_cellularity When a proposed solution of purity and ploidy is tested,
#' Sequenza is run on an interval of purity centered at the proposed value, of half 
#' width \code{delta_cellularity}.
#' @param delta_ploidy When a proposed solution of purity and ploidy is tested,
#' Sequenza is run on an interval of ploidy centered at the proposed value, of half 
#' width \code{delta_ploidy}.
#' @param verbose If TRUE (default FALSE) at each Sequenza fitting step, print the updated
#' list of purity and ploidy pairs to test.
#' @param ... Optional parameters passed to the \code{analyse_peaks} function
#' by CNAqc. Tune these to change error tolerance or karyotypes to use for QC.
#' @return A \code{tibble} containing for each run its number
#' \code{run}, the solution values of \code{purity} and \code{ploidy}, the corresponding
#' CNAqc \code{score} and quality control status \code{QC}, a list \code{sequenza} of the
#' input and output of the Sequenza fitting procedure, and the object of class \code{CNAqc}
#' created by the \code{init()} and \code{analyze_peaks()} functions.
#'
#' @examples
#' \dontrun{
#'
#' # Install Sequenza package from https://cran.r-project.org/web/packages/sequenza/vignettes/sequenza.html
#' install.packages("sequenza")
#'
#' # Outside R
#' pip install sequenza-utils
#'
#' # Process a FASTA file to produce a GC Wiggle track file
#' sequenza−utils gc_wiggle −w 50 --fasta hg19.fa -o hg19.gc50Base.wig.gz
#'
#' # Process BAM and Wiggle files to produce a seqz file
#' sequenza−utils bam2seqz -n normal.bam -t tumor.bam --fasta hg19.fa \
#' -gc hg19.gc50Base.wig.gz -o out.seqz.gz
#'
#' # Post-process by binning the original seqz file:
#' sequenza−utils seqz_binning --seqz out.seqz.gz -w 50 -o out small.seqz.gz
#'
#' # Run this pipeline
#' Sequenza_CNAqc(
#'    sample_id = 'tumour',
#'    seqz_file = 'small.seqz.gz', # Binned file
#'    mutations = dataset$mutations, # If using an external set of mutations
#'    sex = "F", # If female
#'    cellularity = c(0.05, 1),
#'    ploidy = c(1.8, 5.4),
#'    reference = 'GRCh38',
#'    normalization.method = 'median',
#'    window = 1e5,
#'    gamma = 280,
#'    kmin = 300,
#'    min.reads.baf = 50,
#'    min.reads = 50,
#'    min.reads.normal = 15,
#'    max.mut.types = 1,
#'    delta_cellularity = 0.05,
#'    delta_ploidy = 0.25,
#'    verbose = TRUE)
#' )
#' }
#'
#' @export

Sequenza_CNAqc = function(sample_id,
                          seqz_file,
                          mutations = NULL,
                          sex,
                          cellularity = c(0.05, 1),
                          ploidy = c(1.8, 5.4),
                          reference = 'GRCh38',
                          normalization.method = 'median',
                          window = 1e5,
                          gamma = 280,
                          kmin = 300,
                          min.reads.baf = 50,
                          min.reads = 50,
                          min.reads.normal = 15,
                          max.mut.types = 1,
                          delta_cellularity = 0.05,
                          delta_ploidy = 0.25,
                          verbose = FALSE)
{
  
  # Auxiliary function: Sequenza check input parameters
  Sequenza_check_inputs = function(sample_id,
                                   seqzFile,
                                   sex,
                                   cellularity,
                                   ploidy,
                                   reference)
  {
    if (!is.character(sample_id))
      cli::cli_abort("Unrecogniseable sample id, will not proceed.")
    
    cli::cli_alert("Sample id: {.field {sample_id}}")
    
    if (!file.exists(seqzFile))
      cli::cli_abort("File {.field {seqzFile}} does not exist, will not proceed.")
    
    sz = format(object.size(seqzFile), units = "auto")
    cli::cli_alert("Sequenza seqz file: {.value \"{seqzFile}\"} [{.field {sz}}]")
    
    if (!(sex %in% c("M", "F")))
      cli::cli_abort("Sex must be either female (F) or male (M).")
    
    sx = dplyr::case_when(sex == "F" ~ "Female (F)",
                          sex == "M" ~ "Male (M)",)
    cli::cli_alert("Sample sex: {.field {sx}}")
    
    if ((cellularity %>% length() != 2) |
        !(cellularity[1] %>% is.numeric()) |
        !(cellularity[2] %>% is.numeric()))
      cli::cli_abort("Cellularity needs to be a numeric vector with two values.")
    
    cli::cli_alert(
      "Testing cellularity (%) in [{.field {cellularity[1]}}, {.field {cellularity[2]}}]"
    )
    
    if ((ploidy %>% length() != 2) |
        !(ploidy[1] %>% is.numeric()) |
        !(ploidy[2] %>% is.numeric()))
      cli::cli_abort("Ploidy needs to be a numeric vector with two values.")
    
    cli::cli_alert("Testing ploidy (real-valued) in [{.field {ploidy[1]}}, {.field {ploidy[2]}}]")
    
    if (!(reference %>% is.character()) |
        !(reference %in% c("GRCh38", "GRCh37", "hg19")))
      cli::cli_abort("Reference must be any of GRCh38 or hg19/GRCh37.")
    
    cli::cli_alert("Reference genome: {.field {reference}}")
  }
  
  
  #### Check for Sequenza library ####
  if (!require(sequenza)) {
    cli::cli_abort("The R Sequenza package is not installed, will not proceed.")
  }
  
  cli::cli_h1("Sequenza CNA calling wrapper with CNAqc")
  
  Sequenza_check_inputs(
    sample_id = sample_id,
    seqzFile = seqz_file,
    sex = sex,
    cellularity = cellularity,
    ploidy = ploidy,
    reference = reference
  )
  
  #### If male use chrY, if female don't ####
  is_female <-  sex == "F"
  chromosomes <- paste0(c(1:22, "X"))
  # chromosomes = 1:22
  
  if (!as.logical(is_female)) {
    chromosomes <- c(chromosomes, 'Y')
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
    max.mut.types = max.mut.types
    )
  
  # cli::cli_h2("Extraction of required information [{crayon::blue('sequenza.extract')}]")
  cat("\n")
  cli::cli_process_start("Seqz pre-processing [{crayon::blue('sequenza.extract')}]")
  cat("\n")
  
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
  
  run_index = 0
  L_cache = L = NULL
  
  repeat {
    
    if(verbose)
    {
      cli::cli_alert("List of parameters to test")
      print(L)
    }
      
    run_index = run_index + 1
    
    # Handle a new run - special case for the first run
    if (run_index > 1)
    {
      L_head = L %>% filter(row_number() == 1)
      L = L %>% filter(row_number() > 1)
      
      # First run - parameters determined by the user
      cellularity = L_head$cellularity + c(-delta_cellularity, delta_cellularity)
      ploidy = L_head$ploidy + c(-delta_ploidy, delta_ploidy)
    }
    
    # New run - parameters determined by the pipeline
    L_cache_new = sequenza_fit_runner(seqzExt,
                                      mutations,
                                      run = run_index,
                                      is_female,
                                      cellularity,
                                      ploidy,
                                      sample_id,
                                      out_dir = paste0("run_", run_index),
                                      reference = reference)
    
    L_cache = L_cache %>% bind_rows(L_cache_new)
    
    # Obtain proposals
    new_proposals = get_proposals(sequenza = L_cache$sequenza[[nrow(L_cache)]],
                                  cnaqc = L_cache$cnaqc[[nrow(L_cache)]])
    
    # Check them against the L_cache, create a list to go for
    new_proposals = filter_proposals(
      L_cache = L_cache,
      new_proposals = new_proposals,
      delta_c = delta_cellularity,
      delta_p = delta_ploidy
    )
    
    
    if(verbose)
    {
      cli::cli_alert("New proposal")
      print(new_proposals)
      
      cli::cli_alert("Cached evaluations")
      print(L_cache)
      
    }
    
    L = L %>% bind_rows(new_proposals)
    
    if (L %>% nrow() == 0)
      break
  }
  
  # End of pipeline
  best_fit = which.min(L_cache$score %>% abs())
  qc_status = L_cache$QC[best_fit]
  
  cli::cli_h2("Best fit with score {.field {L_cache$score[best_fit]}} and QC {.field {qc_status}}")
  
  if(qc_status == "FAIL")
    cli::cli_alert_danger("")
  else
    cli::cli_alert_success("")
  
  R.utils::createLink(
    link = 'final',
    target = L_cache$run[best_fit],
    skip = FALSE,
    overwrite = TRUE
  )
  
  # Plot all CNAqc results
  pdf("./cnaqc_reports.pdf", onefile = TRUE, width = 18, height = 10)  
  lapply(L_cache$cnaqc %>% seq_along, function(x){
    ggpubr::ggarrange(
      plot_segments(L_cache$cnaqc[[x]], highlight = c('2:0', '2:1', '2:2', '1:1', '1:0'))+
        ggtitle(paste("Run", x)),
      plot_peaks_analysis(L_cache$cnaqc[[x]]),
      ncol = 1) %>% print()
  })
  dev.off()

  # Save cumulative pipeline results
  saveRDS(L_cache, file = "pipeline.rds")
  
  # Return summary table
  return(L_cache)
}

## Auxiliary function: Sequenza extract and fit data
sequenza_fit_runner = function(seqzExt,
                               mutations,
                               run,
                               is_female,
                               cellularity,
                               ploidy,
                               sample_id,
                               out_dir,
                               reference)
{
  # Auxiliary function: Sequenza outputs parser
  parse_Sequenza_CNAs = function(out, x)
  {
    # Files required
    segments_file = paste0(out, '/', x, "_segments.txt")
    mutations_file = paste0(out, '/', x, "_mutations.txt")
    purity_file = paste0(out, '/', x, "_confints_CP.txt")
    alternatives_file = paste0(out, '/', x, "_alternative_solutions.txt")
    
    if (!file.exists(segments_file))
      stop(paste0("Missing required segments file: ", segments_file))
    
    if (!file.exists(mutations_file))
      stop(paste0("Missing required mutations file: ", mutations_file))
    
    if (!file.exists(purity_file))
      stop(paste0("Missing required solutions file: ", purity_file))
    
    if (!file.exists(alternatives_file))
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
      dplyr::rename(chr = chromosome,
                    from = position,) %>%
      dplyr::rowwise() %>%
      tidyr::separate(
        mutation,
        sep = '>',
        into = c('ref', 'alt'),
        remove = FALSE
      ) %>%
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
  
  
  if (dir.exists(out_dir))
    cli::cli_alert_warning("Directory {.field {out_dir}} already exists")
  
  dir.create(out_dir)
  
  chr.fit = seqzExt$chromosomes %>% unique()
  # Select only the more reliable autosomes for fitting cellularity and ploidy parameters
  chr.fit = chr.fit[!(chr.fit %in% c('X', 'Y', 'chrX', 'chrY'))]
  
  cat("\n")
  cli::cli_process_start(
    "Sequenza fit via CNAqc [{crayon::blue('sequenza.fit')}] - cellularity {.value {cellularity}}, ploidy {.value {ploidy}}"
  )
  cat("\n")
  
  # Run of sequenza.fit
  if (cellularity[1] < 0)
    cellularity[1] = 0.01
  if (cellularity[1] > 1)
    cellularity[1] = 1
  
  if (cellularity[2] < 0)
    cellularity[2] = 0.01
  if (cellularity[2] > 1)
    cellularity[2] = 1
  
  if (ploidy[1] < 0)
    ploidy[1] = 0.01
  if (ploidy[2] < 0)
    ploidy[2] = 0.01
  
  if(run == 1){
    nbins_cellularity = 100
    nbins_ploidy = 50
  }
  else{
    nbins_cellularity = 10
    nbins_ploidy = 10
  }
  paraSpace <- sequenza::sequenza.fit(
    sequenza.extract = seqzExt,
    cellularity = seq(cellularity[1], cellularity[2], length.out = nbins_cellularity),
    ploidy = seq(ploidy[1], ploidy[2], length.out = nbins_ploidy),
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
  
  cat("\n")
  cli::cli_process_done()
  
  # Parse outputs
  cat("\n")
  cli::cli_process_start("Quality control with CNAqc: {.field {out_dir}}")
  cat("\n")
  
  fits = parse_Sequenza_CNAs(out_dir, x = sample_id)
  
  # Perform QC, and save RDS
  cnaqc_obj = tryCatch({
    if(is.null(mutations)){
      cat("\n")
      cli::cli_alert("Using Sequenza mutation calls")
      cat("\n")
      mutations = fits$mutations %>%
        rename(VAF = "F", DP = "good.reads") %>%
        mutate(NV = as.integer(VAF * DP))
    }
    else{
      cat("\n")
      cli::cli_alert("Using the input set {.field {mutations}} as mutation calls")
      cat("\n")
    }
    CNAqc::init(
      snvs =  mutations,
      cna = fits$segments,
      purity = fits$purity,
      ref = reference
    ) %>%
      CNAqc::analyze_peaks()
  },
  error = function(e)
  {
    print(e)
    cli::cli_abort("CNAqc generated errors and cannot run on this sample. Aborting the pipeline.")
  })
  
  purity = fits$purity[1]
  ploidy = fits$ploidy[1]
  
  # show to screen
  print(cnaqc_obj)
  
  cat("\n")
  cli::cli_process_done()
  
  
  # L_cache = append(L_cache, list(c(purity,ploidy)))
  tibble(
    run = out_dir,
    purity = purity,
    ploidy = ploidy,
    score = cnaqc_obj$peaks_analysis$score,
    QC = cnaqc_obj$peaks_analysis$QC,
    sequenza = list(fits),
    cnaqc = list(cnaqc_obj)
  )
}

get_proposals = function(sequenza, cnaqc)
{
  sqz = sequenza$alternative_solutions %>% select(-SLPP)
  
  cnq = tibble(
    cellularity = cnaqc$purity + cnaqc$peaks_analysis$score,
    ploidy = sequenza$ploidy
  )
  
  bind_rows(sqz,
            cnq) %>%
    dplyr::filter(cellularity <= 1 & cellularity > 0)
}

filter_proposals = function(L_cache,
                            new_proposals,
                            delta_cellularity,
                            delta_ploidy)
{
  bitmask = c()
  
  for (i in 1:nrow(new_proposals))
  {
    cellularity_cache = abs(L_cache$purity - new_proposals$cellularity[i]) > delta_cellularity
    ploidy_cache = abs(L_cache$ploidy - new_proposals$ploidy[i]) > delta_ploidy
    
    bitmask = c(bitmask, all(cellularity_cache | ploidy_cache))
  }
  
  new_proposals %>%
    dplyr::filter(bitmask)
}
