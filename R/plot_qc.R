#' Plot a summary of QC results.
#'
#'  @description Results from \code{analyze_peaks} can be visualised with this
#'  function. Compared to individual karyotypes fits available with function \code{plot_peaks_analysis},
#'  this function reports sumary statistics for each karyotype, and the overall score.
#'
#' @param x An object of class \code{cnaqc}, where function \code{analyze_peaks} has
#' been computed.
#' @param digits Number of digits used to round scores.
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna,example_dataset_CNAqc$purity)
#'
#' x = analyze_peaks(x)
#' plot_qc(x)
plot_qc = function(x, digits = 4)
{
  stopifnot(inherits(x, "cnaqc"))

  with_peaks = all(!is.null(x$peaks_analysis))
  if(!with_peaks){
    stop("Missing peaks, see ?peaks_analysis to run peaks analysis.")
  }

  with_ccf = all(!is.null(x$CCF_estimates))
  if(!with_ccf){
    stop("Missing CCF, see ?compute_CCF to estimate CCF values.")
  }

  all_karyptypes = x$peaks_analysis$plots %>% names

  button = function(color, radio)
  {
    ggplot(data = data.frame(
      x = 0,
      y = 0,
      label = radio,
      stringsAsFactors = FALSE
    ),
    aes(x = x , y = y, label = label)) +
      theme_void() +
      theme(plot.background = element_rect(fill = color)) +
      geom_text()
  }

  # peak detection
  analysed_karyotypes = x$peaks_analysis$matches$karyotype %>% unique
  qc_karyotypes = x$peaks_analysis$matches %>%
    dplyr::distinct(karyotype, QC)

  # create remote
  radios_peaks = lapply(all_karyptypes,
         function(k) {
           if (k %in% analysed_karyotypes)
           {
             QC = qc_karyotypes %>% dplyr::filter(karyotype == k) %>% dplyr::pull(QC)
             if(QC == "PASS") return(button('forestgreen', k))
             else return(button('indianred3', k))
           }
           else
             return(button('gainsboro', k))
         })

  radios_peaks = append(list(ggpubr::text_grob("Peaks", face = 'bold')), radios_peaks)
  radios_peaks = ggpubr::ggarrange(
    plotlist = radios_peaks, ncol = length(radios_peaks), nrow = 1)

  # CCF
  analysed_karyotypes = x$CCF_estimates %>% names

  # create remote
  radios_ccf = lapply(all_karyptypes,
                        function(k) {
                          if (k %in% analysed_karyotypes)
                          {
                            QC = x$CCF_estimates[[k]]$QC_table$QC
                            if(QC == "PASS") return(button('forestgreen', k))
                            else return(button('indianred3', k))
                          }
                          else
                            return(button('gainsboro', k))
                        })

  radios_ccf = append(list(ggpubr::text_grob("CCF", face = 'bold')), radios_ccf)
  radios_ccf = ggpubr::ggarrange(
    plotlist = radios_ccf, ncol = length(radios_ccf), nrow = 1)


  remotes =
    ggpubr::ggarrange(radios_peaks, radios_ccf, nrow = 2, ncol = 1)

  return(remotes)
}
