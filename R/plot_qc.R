#' Plot a summary of QC results.
#'
#' @description Results from \code{analyze_peaks} and \code{compute_CCF} can be
#' visualised with this function. Compared to individual karyotypes fits available
#' with function \code{plot_peaks_analysis}, for instance, this function reports
#' summary pass/fail statistics for each analysis.
#'
#' @param x A CNAqc object.
#'
#' @return A \code{ggplot2} plot
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(mutations = example_dataset_CNAqc$mutations, cna = example_dataset_CNAqc$cna, purity = example_dataset_CNAqc$purity)
#'
#' x = analyze_peaks(x)
#' x = compute_CCF(x)
#'
#' plot_qc(x)
plot_qc = function(x)
{
  # stopifnot(inherits(x, "cnaqc"))
  #
  # with_peaks = all(!is.null(x$peaks_analysis))
  # if(!with_peaks){
  #   stop("Missing peaks, see ?peaks_analysis to run peaks analysis.")
  # }
  #
  # with_ccf = all(!is.null(x$CCF_estimates))
  # if(!with_ccf){
  #   stop("Missing CCF, see ?compute_CCF to estimate CCF values.")
  # }
  #
  # all_karyptypes = x$peaks_analysis$plots %>% names
  #
  # button = function(color, radio)
  # {
  #   ggplot(data = data.frame(
  #     x = 0,
  #     y = 0,
  #     label = radio,
  #     stringsAsFactors = FALSE
  #   ),
  #   aes(x = x , y = y, label = label)) +
  #     theme_void() +
  #     theme(plot.background = element_rect(fill = color)) +
  #     geom_text()
  # }
  #
  # # peak detection
  # analysed_karyotypes = x$peaks_analysis$matches$karyotype %>% unique
  # qc_karyotypes = x$peaks_analysis$matches %>%
  #   dplyr::distinct(karyotype, QC)
  #
  # # create remote
  # radios_peaks = lapply(all_karyptypes,
  #        function(k) {
  #          if (k %in% analysed_karyotypes)
  #          {
  #            QC = qc_karyotypes %>% dplyr::filter(karyotype == k) %>% dplyr::pull(QC)
  #            if(QC == "PASS") return(button('forestgreen', k))
  #            else return(button('indianred3', k))
  #          }
  #          else
  #            return(button('gainsboro', k))
  #        })
  #
  # radios_peaks = append(list(ggpubr::text_grob("Peaks", face = 'bold')), radios_peaks)
  # radios_peaks = ggpubr::ggarrange(
  #   plotlist = radios_peaks, ncol = length(radios_peaks), nrow = 1)
  #
  # # CCF
  # analysed_karyotypes = x$CCF_estimates %>% names
  #
  # # create remote
  # radios_ccf = lapply(all_karyptypes,
  #                       function(k) {
  #                         if (k %in% analysed_karyotypes)
  #                         {
  #                           QC = x$CCF_estimates[[k]]$QC_table$QC
  #                           if(QC == "PASS") return(button('forestgreen', k))
  #                           else return(button('indianred3', k))
  #                         }
  #                         else
  #                           return(button('gainsboro', k))
  #                       })
  #
  # radios_ccf = append(list(ggpubr::text_grob("CCF", face = 'bold')), radios_ccf)
  # radios_ccf = ggpubr::ggarrange(
  #   plotlist = radios_ccf, ncol = length(radios_ccf), nrow = 1)
  #
  #
  # remotes =
  #   ggpubr::ggarrange(radios_peaks, radios_ccf, nrow = 2, ncol = 1)
  #
  # return(remotes)

  # stopifnot(inherits(x, "cnaqc"))
  #
  # with_peaks = all(!is.null(x$peaks_analysis))
  # if(!with_peaks){
  #   stop("Missing peaks, see ?peaks_analysis to run peaks analysis.")
  # }
  #
  # with_ccf = all(!is.null(x$CCF_estimates))
  # if(!with_ccf){
  #   stop("Missing CCF, see ?compute_CCF to estimate CCF values.")
  # }
  #
  # all_karyptypes = x$peaks_analysis$plots %>% names
  #
  # peaks_QC = x$peaks_analysis$matches %>%
  #   dplyr::select(karyotype, QC) %>%
  #   dplyr::full_join(data.frame(karyotype = all_karyptypes, stringsAsFactors = F), by = 'karyotype') %>%
  #   dplyr::distinct(karyotype, QC, .keep_all = T) %>%
  #   dplyr::arrange(karyotype) %>%
  #   dplyr::mutate(
  #     value = 1,
  #     lab.ypos = cumsum(value) - 0.5 * value,
  #     QC = paste(QC),
  #     label = karyotype,
  #     type = 'Peaks')
  #
  # CCF_QC = Reduce(dplyr::bind_rows, lapply(x$CCF_estimates, function(x) x$QC_table)) %>%
  #   dplyr::select(karyotype, QC) %>%
  #   dplyr::full_join(data.frame(karyotype = all_karyptypes, stringsAsFactors = F), by = 'karyotype') %>%
  #   dplyr::arrange(karyotype) %>%
  #   dplyr::mutate(
  #     value = 1,
  #     lab.ypos = cumsum(value) - 0.5 * value,
  #     QC = paste(QC),
  #     label = karyotype,
  #     type = 'CCF')
  #
  # QC_table = dplyr::bind_rows(peaks_QC, CCF_QC)
  # QC_table$karyotype = factor(QC_table$karyotype, all_karyptypes)
  # QC_table$type = factor(QC_table$type, levels = c('Peaks', 'CCF'))
  #
  # QC_table %>%
  #   ggplot(aes(x = 1, y = value, fill = QC)) +
  #   facet_wrap(~type) +
  #   geom_bar(width = 1,
  #            stat = "identity",
  #            color = "white") +
  #   coord_polar("y", start = 0) +
  #   geom_label(aes(y = lab.ypos, label = karyotype),
  #              color = "black",
  #              fill = 'white') +
  #   scale_fill_manual(values = c(
  #     `PASS` = 'forestgreen',
  #     `FAIL` = 'indianred3',
  #     `NA` = 'gainsboro'
  #   )) +
  #   CNAqc:::my_ggplot_theme() +
  #   theme(
  #     axis.text.x = element_blank(), axis.ticks.x = element_blank(),
  #     axis.text.y = element_blank(), axis.ticks.y = element_blank()
  #   ) +
  #   labs(x = NULL, y = NULL, title = "CNAqc summary QC")

  stopifnot(inherits(x, "cnaqc"))

  xqc = compute_QC_table(x)

  QC_table = xqc$QC_table
  pPASS = xqc$percentage_PASS
  NA_tests = xqc$NA_tests

  top_table = QC_table %>%
    ggplot2::ggplot(aes(x = karyotype, y = type, fill = paste(QC))) +
    ggplot2::facet_wrap( ~ paste0("QC (% of PASS is ", pPASS, '%, - NA tests are ', NA_tests, '%)')) +
    ggplot2::geom_tile(aes(width = .8, height = .8)) +
    ggplot2::scale_fill_manual(values = c(
      `PASS` = 'forestgreen',
      `FAIL` = 'indianred3',
      `NA` = 'gainsboro'
    )) +
    CNAqc:::my_ggplot_theme() +
    ggplot2::labs(x = NULL, y = NULL, title = "Simple clonal CNAs ") +
    ggplot2::guides(fill = ggplot2::guide_legend('QC (NA not available)'))

  # Other peaks
  secondary_table = ggplot2::ggplot() + ggplot2::labs(title = "Complex clonal CNAs")
  if (!is.null(x$peaks_analysis$general))
    secondary_table = secondary_table +
    ggplot2::geom_tile(
      data = x$peaks_analysis$general$expected_peaks,
      ggplot2::aes(
        y = karyotype,
        x = multiplicity,
        fill = matched,
        width = .8,
        height = .8
      )
    ) +
    CNAqc:::my_ggplot_theme() +
    ggplot2::scale_fill_manual(values = c(`FALSE` = 'indianred3', `TRUE` = 'forestgreen')) +
    ggplot2::facet_wrap( ~ "General peaks")

  # Subclonal peaks
  third_table = ggplot() + labs(title = "Subclonal CNAs")
  if (!is.null(x$peaks_analysis$subclonal))
  {
    expected_peaks = x$peaks_analysis$subclonal$expected_peaks

    third_table = ggplot2::ggplot(expected_peaks %>% group_by(segment_id, model) %>% mutate(peak = row_number())) +
      ggplot2::geom_tile(
        ggplot2::aes(
          y = segment_id,
          x = peak,
          fill = matched,
          color = role
        ),
        width = .8,
        size = 1,
        height = .8
      ) +
      ggplot2::scale_fill_manual(values = c(
        `TRUE` = 'forestgreen',
        `FALSE` = 'indianred',
        `NA` = 'gray'
      )) +
      ggplot2::facet_wrap( ~ model) +
      CNAqc:::my_ggplot_theme() +
      ggplot2::scale_y_discrete(limits = expected_peaks$segment_id %>% unique %>% gtools::mixedsort(decreasing = TRUE)) +
      # scale_x_discrete(limits = c(1:2) %>% paste) +
      ggplot2::guides(
        fill = ggplot2::guide_legend(ncol = 1),
        color = ggplot2::guide_legend(override.aes = ggplot2::aes(fill = NA), ncol = 1)
      ) +
      ggplot2::labs(y = 'Subclonal CNA segment') +
      ggplot2::scale_color_manual(values = c(`private` = 'gray', `shared` = 'black')) +
      ggplot2::coord_flip() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) +
      ggplot2::labs(title = 'Subclonal CNAs')
  }

  p1 = cowplot::plot_grid(
    top_table,
    secondary_table,
    axis = 'tb',
    align = 'h',
    ncol = 2
  )
  cowplot::plot_grid(p1, third_table, ncol = 1)

  # cowplot::plot_grid(top_table, secondary_table, third_table, axis = 'tb', align = 'h', ncol = 3)
}


compute_QC_table = function(x)
{
  all_karyptypes = c("1:0", "1:1", "2:1", "2:0", '2:2')

  # Tables for peaks
  if (all(is.null(x$peaks_analysis)))
  {
    peaks_QC = data.frame(
      karyotype = all_karyptypes,
      QC = NA,
      type = 'Peaks',
      stringsAsFactors = FALSE
    ) %>% as_tibble()
  }
  else
  {
    peaks_QC = x$peaks_analysis$matches %>%
      dplyr::select(karyotype, QC) %>%
      dplyr::full_join(data.frame(karyotype = all_karyptypes, stringsAsFactors = F),
                       by = 'karyotype') %>%
      dplyr::distinct(karyotype, QC, .keep_all = T) %>%
      dplyr::arrange(karyotype) %>%
      dplyr::mutate(# value = 1,
        # lab.ypos = cumsum(value) - 0.5 * value,
        QC = paste(QC),
        # label = karyotype,
        type = 'Peaks')
  }

  # Table for CCF
  if (all(is.null(x$CCF_estimates)))
  {
    CCF_QC = data.frame(
      karyotype = all_karyptypes,
      QC = NA,
      type = 'CCF',
      stringsAsFactors = FALSE
    ) %>% as_tibble()
  }
  else
  {
    CCF_QC = Reduce(dplyr::bind_rows,
                    lapply(x$CCF_estimates, function(x)
                      x$QC_table)) %>%
      dplyr::select(karyotype, QC) %>%
      dplyr::full_join(data.frame(karyotype = all_karyptypes, stringsAsFactors = F),
                       by = 'karyotype') %>%
      dplyr::arrange(karyotype) %>%
      dplyr::mutate(# value = 1,
        # lab.ypos = cumsum(value) - 0.5 * value,
        QC = paste(QC),
        # label = karyotype,
        type = 'CCF')
  }

  # Bind both tables
  QC_table = dplyr::bind_rows(peaks_QC, CCF_QC) %>%
    dplyr::mutate(QC = ifelse(!is.na(QC) & QC == "NA", NA, QC))

  # QC_table$karyotype = factor(QC_table$karyotype, all_karyptypes)
  # QC_table$type = factor(QC_table$type, levels = c('Peaks', 'CCF'))
  # QC_table = QC_table %>%
  #   dplyr::mutate(
  #     QC = ifelse(!is.na(QC) & QC == "NA", NA, QC)
  #   )

  # Percentage of PASS cases, out of not-NA ones
  pPASS =
    sum(QC_table$QC == "PASS", na.rm = T) / (sum(QC_table$QC == "PASS", na.rm = T) + sum(QC_table$QC == "FAIL", na.rm = T)) * 100

  # Total number of NA tests
  NA_tests = sum(is.na(QC_table$QC)) / nrow(QC_table) * 100

  return(list(
    QC_table = QC_table,
    percentage_PASS = pPASS,
    NA_tests = NA_tests
  ))
}
