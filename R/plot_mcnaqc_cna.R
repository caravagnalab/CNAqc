#' Plot multisample cna data
#'
#' @param x A m_cnaqc object 
#' @param which either "original" or "shared"
#' @param chromosomes 
#'
#' @return a plot
#' @export
#'
#' @examples 
#' \dontrun{plot_segments_multisample(x, which = "shared")}
plot_segments_multisample <- function(x, 
                                      which, 
                                      chromosomes = paste0('chr', c(1:22, 'X', 'Y'))) {
  # add here a sanitasier
  
  Ln = get_sample_name(x)
  if (is.null(Ln)) {
    Ln = paste0("Sample ", 1:length(x))
    names(x) = Ln
    
    cli::cli_alert_warning("The input list is un-named, using default naming scheme Sample*")
  }
  
  KARYO_colors = CNAqc:::get_karyotypes_colors(NULL)
  
  cnaqc_list = get_sample(x, sample = Ln, which_obj = which)
  
  # Extract calls, and flatten them for plotting --> on mcnaqc obj
  calls = lapply(cnaqc_list,
                 function(s)
                 {
                   W = CNA(s) %>%
                     filter(chr %in% chromosomes) %>% 
                     mutate(
                       label = paste(Major, minor, sep = ':'),
                       CN = minor + Major,
                       sample_id = get_sample_name(s)
                     ) %>%
                     select(chr, from, to, label, CN, sample_id)
                   
                   CNAqc:::relative_to_absolute_coordinates(s, W)
                 })
  
  calls_flat =
    suppressWarnings(Reduce(function(x, y)
      full_join(
        x, y, by = c("chr", "from", "to", "label", "CN", "sample_id")
      ),
      calls) %>%
        mutate(label = ifelse(
          label %in% names(KARYO_colors), label, 'other'
        )))
  
  KARYO_colors = c(KARYO_colors, `other` = 'gray')
  
  # chromosomes = calls_flat$chr %>% unique
  
  # Reference genome
  ref_gen = cnaqc_list[[1]]$reference_genome
  ref_coords =cnaqc_list[[1]]$genomic_coordinates
  reference_genome = CNAqc:::get_reference(ref_gen, data = ref_coords) %>% 
    dplyr::filter(chr %in% chromosomes)
  low = min(reference_genome$from)
  upp = max(reference_genome$to)
  
  bl_genome = suppressMessages(
    blank_genome_multisample(x,
      ref = ref_gen,
      reference_coordinates = reference_genome,
      chromosomes = chromosomes,
      label_chr = NA
    ) +
      theme(axis.title.y = element_text(margin = margin(r = 20)), 
            axis.title.x = element_text(margin = margin(t = 15)))
  )
  
  cna_multisample = bl_genome +
    geom_linerange(data = calls_flat, 
                   aes(xmin = from, 
                       xmax = to, 
                       colour = label), size = 5) +
    ggplot2::scale_color_manual(values = KARYO_colors) +
    ggplot2::ggtitle("Comparative CNA")
    
  
  return(cna_multisample)
}
