#' Print for class \code{'cnaqc'}.
#'
#' @param x An obj of class \code{'cnaqc'}.
#' @param ... Default S3 method parameter.
#'
#' @return Nothing.
#'
#' @exportS3Method print cnaqc
#' @export print.cnaqc
#'
#' @import pio
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna,example_dataset_CNAqc$purity)
#'
#' print(x)
print.cnaqc = function(x, ...)
{
  stopifnot(inherits(x, "cnaqc"))

  cli::cli_rule(
    paste(
      crayon::bgYellow(crayon::black("[ CNAqc ] ")),
      'n = {.value {x$n_snvs}} mutations in {.field {x$n_cna}} segments ({.value {x$n_cna_clonal}} clonal + {.value {x$n_cna_sbclonal}} subclonal). Genome reference: {.field {x$reference_genome}}.'
    )
  )

  # cli::cli_alert_info(paste0(" CNA segments: ", x$n_cna_clonal, " clonal, ", x$n_cna_sbclonal, " subclonal."))

  # cli::cli_alert_info(paste0("Mutation mapping (head): ", paste0(head(x$n_karyotype), ' (',
  # names(head(x$n_karyotype)), ')', collapse = '; ')))

  # cli::cli_alert_info(paste0("Mutation mapping (up to top 5): "))
  cat('\n')
  bar_print_console(x)
  cat('\n')

  # Available analyses
  with_peaks = all(!is.null(x$peaks_analysis))
  with_CCF = all(!is.null(x$CCF_estimates))
  with_smoothing = all(!is.null(x$before_smoothing))
  with_arm_frag = all(!is.null(x$arm_fragmentation))
  with_wg_frag = all(!is.null(x$wg_fragmentation))
  with_drivers = all(c("gene", "is_driver") %in% colnames(x$snvs))

  cli::cli_alert_info(paste0(
    "Sample Purity: ",
    paste0(x$purity  * 100, '% ~ Ploidy: ', x$ploidy, '.')
  ))

  if (with_drivers |
      with_peaks |
      with_CCF | with_smoothing | with_arm_frag | with_wg_frag)
    cat('\n')

  # if(with_drivers)
  # {
  #   nd = x$snvs %>% dplyr::filter(is_driver) %>% nrow()
  #   cli::cli_alert_info("Mutations annotated have {.value {nd}} drivers.")
  # }

  ppass = function()
    "{crayon::bgGreen(crayon::black(\" PASS \"))}"
  pfail = function()
    "{crayon::bgRed(crayon::white(\" FAIL \"))}"

  if (with_peaks)
  {
    prop = x$peaks_analysis$matches %>%
      dplyr::group_by(QC) %>%
      dplyr::summarise(prop = sum(weight), .groups = 'drop') %>%
      dplyr::arrange(desc(prop)) %>%
      dplyr::filter(dplyr::row_number() == 1) %>%
      dplyr::pull(prop)
    prop = round(prop * 100, digits = 0)

    if (x$peaks_analysis$QC == "PASS")
      cli::cli_h3(
        paste0(
          ppass(),
          " Peaks QC {crayon::bold(x$peaks_analysis$matching_strategy)}: % {crayon::green(paste0(prop))} with {crayon::green(paste('q =', x$peaks_analysis$score))}"
        )
      )
    else
      cli::cli_h3(
        paste0(
          pfail(),
          " Peaks QC {crayon::bold(x$peaks_analysis$matching_strategy)}: % {crayon::red(paste0(prop))} with {crayon::red(paste('q =', x$peaks_analysis$score))}"
        )
      )

    xx = x$peaks_analysis$matches %>%
      dplyr::mutate(QC = ifelse(QC == "PASS", ppass(), pfail()))

    for (karyo in xx$karyotype %>% unique)
    {
      xxx = xx %>%
        dplyr::filter(karyotype == karyo)

      qc = paste(xxx$QC, sprintf("%-7s", round(xxx$offset, 3)), collapse = ' ')
      # qc = paste0("[", xx$karyotype[1], "]",  qc)

      n = sprintf("%-5s", x$n_karyotype[karyo])
      p = x$n_karyotype[karyo] / sum(x$n_karyotype[xx$karyotype %>% unique])
      p = format(p * 100, digits = 0)
      p = sprintf("%2s", p)

      cli::cli_alert_info(
        paste0(
          crayon::blue(xxx$karyotype[1]),
          " ~ n = {n} ({p}%) {clisymbols::symbol$arrow_right} ",
          qc,
          ""
        )
      )

    }
  }

  if (with_CCF)
  {
    cat("\n")
    cli::cli_alert_success(
      "Cancer Cell Fraction (CCF) data available for karyotypes:{.value {names(x$CCF_estimates)}}."
    )

    lapply(x$CCF_estimates,
           function(l) {
             l$QC_table
           }) %>%
      Reduce(f = bind_rows) %>%
      apply(1, function(z) {
        if (z['QC'] == "PASS")
          cli::cli_alert_success(
            "{
            crayon::bgGreen(
              crayon::white(
                \" PASS \"))} CCF via {.value {crayon::bold(z['method'])}}."
          )
        else
          cli::cli_alert_success(
            "{crayon::bgRed(crayon::white(\" FAIL \"))} CCF via {.value {crayon::bold(z['method'])}}."
          )
      })
  }



  if (with_smoothing)
    cli::cli_alert_success(
      "These segments are smoothed; before smoothing there were {.value {x$before_smoothing$n_cna}} segments."
    )

  if (with_arm_frag)
    cli::cli_alert_success(
      "Arm-level fragmentation analysis: {.value {sum(x$arm_fragmentation$table$significant)}} segments overfragmented."
    )

  if (with_wg_frag)
  {
    cond = x$wg_fragmentation$is_overfragmented
    p = round(x$wg_fragmentation$pvalue, 5)
    cli::cli_alert_success(
      "Whole-genome fragmentation analysis: p = {.value {p}}: {.value {ifelse(cond, crayon::red('overfragmented'), crayon::green('not overfragmented'))}}."
    )
  }

  ppass = function()
    "{crayon::bgGreen(crayon::black(\" PASS \"))}"
  pfail = function()
    "{crayon::bgRed(crayon::white(\" FAIL \"))}"

}

bar_print_console = function(x, top = length(x$n_karyotype)) {
  e = x$n_karyotype %>% sort(decreasing = TRUE)
  l =  round(x$l_karyotype / 10 ^ 6, digits = 0)


  width = options("width")$width / 3 %>% round
  bars = (e / max(e) * width) %>% round

  bars = bars[1:min(top, length(bars))]

  # Max labels widts
  label_width = sapply(names(e), nchar) %>% max
  entries_width = sapply(e, nchar) %>% max

  p = sapply(names(bars),
             function(b)
             {
               cat(sprintf(' %*s ', label_width, b))
               cat(sprintf(' [n = %*s, L = %*s Mb] ',
                           entries_width,
                           e[b],
                           entries_width,
                           l[b]))

               cat(paste(rep("\u25A0", bars[b]), collapse = ''))

               drv = x$snvs %>% dplyr::filter(karyotype == b)
               if ('is_driver' %in% colnames(x$snvs))
               {
                 drv = drv %>% dplyr::filter(is_driver) %>% dplyr::pull(driver_label)

                 if (length(drv) > 0)
                   cat('  {', crayon::yellow(paste(drv, collapse = ', ')), '}')
               }

               cat("\n")
             })
}


#' Plot for class \code{'cnaqc'}.
#'
#' @description
#'
#' The default plot is the CNA segments in wide format
#'
#' @param x An obj of class \code{'cnaqc'}.
#' @param ... Default S3 method parameter.
#'
#' @return Nothing.
#'
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna,example_dataset_CNAqc$purity)
#'
#' plot.cnaqc(x)
plot.cnaqc = function(x, ...)
{
  stopifnot(inherits(x, "cnaqc"))

  with_CCF = all(!is.null(x$CCF_estimates))

  # Annotate genome wide muts
  top_plots = suppressWarnings(suppressMessages(
    list(
      # plot_gw_counts(x),
      # plot_gw_vaf(x, N = 10000),
      plot_gw_depth(x, N = 10000)
      # plot_gw_ccf(x, N = 10000)
    )
  ))

  # Segment plots
  segments_plot = cowplot::plot_grid(
    plotlist = append(top_plots, list(plot_segments(x))),
    align = 'v',
    nrow = length(top_plots) + 1,
    rel_heights = c(rep(.15, length(top_plots)), .3 * length(top_plots) + 0.2)
  )

  # Data histograms
  top_plots = suppressWarnings(suppressMessages(
    list(
      plot_data_histogram(x, which = 'CCF'),
      plot_data_histogram(x, which = 'VAF'),
      plot_data_histogram(x, which = 'DP'),
      plot_data_histogram(x, which = 'NV')
    )
  ))

  hist_plot = ggpubr::ggarrange(
    plotlist = top_plots,
    nrow = 1,
    ncol = 4,
    common.legend = TRUE,
    legend = 'bottom'
  )

  # Figure
  ggpubr::ggarrange(
    segments_plot,
    hist_plot,
    nrow = 2,
    ncol = 1,
    heights = c(1, .7)
  )



}
