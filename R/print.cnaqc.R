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
#' @examples
#' \dontrun{
#' data('example_PCAWG', package = 'CNAqc')
#'
#' print(example_PCAWG)
#' }
print.cnaqc = function(x, ...)
{
  stopifnot(inherits(x, "cnaqc"))

  cli::cli_rule(
    paste(
      crayon::bgYellow(crayon::black(paste0("[ CNAqc ] ", x$sample))),
      '{.field {x$n_mutations}} mutations in {.field {x$n_cna}} segments ({.field {x$n_cna_clonal}} clonal, {.field {x$n_cna_subclonal}} subclonal). Genome reference: {.field {x$reference_genome}}.'
    )
  )

  # cli::cli_alert_info(paste0(" CNA segments: ", x$n_cna_clonal, " clonal, ", x$n_cna_sbclonal, " subclonal."))

  # cli::cli_alert_info(paste0("Mutation mapping (head): ", paste0(head(x$n_karyotype), ' (',
  # names(head(x$n_karyotype)), ')', collapse = '; ')))

  # cli::cli_alert_info(paste0("Mutation mapping (up to top 5): "))

  cli::cli_h3("Clonal CNAs")
  cat('\n')
  bar_print_console(x, top = 10)
  cat('\n')

  if(x$n_cna_subclonal > 0)
  {
    cli::cli_h3("Subclonal CNAs (showing up to 10 segments)")
    cat('\n')
    bar_print_console_scl(x, top = 10)
    cat('\n')
  }

  # Available analyses
  with_peaks = all(!is.null(x$peaks_analysis))
  with_CCF = all(!is.null(x$CCF_estimates))
  with_smoothing = all(!is.null(x$before_smoothing))
  with_arm_frag = all(!is.null(x$arm_fragmentation))
  with_wg_frag = all(!is.null(x$wg_fragmentation))
  with_drivers = x %>% has_driver_data()
  with_MAF = x %>% has_MAF_annotations()

  cli::cli_alert_info(paste0(
    "Sample Purity: ",
    paste0(x$purity  * 100, '% ~ Ploidy: ', x$ploidy, '.')
  ))

  if (with_drivers |
      with_peaks |
      with_CCF | with_smoothing | with_arm_frag | with_wg_frag | with_MAF)
    cat('\n')

  if(with_drivers)
  {
    nd = x$mutations %>% dplyr::filter(is_driver) %>% nrow()
    cli::cli_alert_info("There are {.value {nd}} annotated driver(s) mapped to clonal CNAs.")

    w_d = x$mutations %>%
      dplyr::filter(is_driver) %>%
      dplyr::select(chr, from, to, ref, alt, DP, NV, VAF, driver_label, is_driver) %>%
      as.data.frame()

    writeLines(paste0("      ",
                      capture.output( w_d %>%
                          print(row.names = F)
                      )))
  }

  if(with_MAF)
  {
    cli::cli_alert("This sample seem to have MAF columns annotated")
  }

  ppass = function()
    "{crayon::bgGreen(crayon::black(\" PASS \"))}"
  pfail = function()
    "{crayon::bgRed(crayon::white(\" FAIL \"))}"


  if (with_peaks)
  {
    if ("matches" %in% names(x$peaks_analysis))
    {
      prop = x$peaks_analysis$matches %>%
        dplyr::group_by(QC) %>%
        dplyr::summarise(prop = sum(weight), .groups = 'drop') %>%
        dplyr::arrange(dplyr::desc(prop)) %>%
        dplyr::filter(dplyr::row_number() == 1) %>%
        dplyr::pull(prop)
      prop = round(prop * 100, digits = 0)

      pur_sc = round(x$peaks_analysis$score, digits = 4)
      pur_ch = paste0(round(x$peaks_analysis$score * 100, digits = 0), '%')

      if (x$peaks_analysis$QC == "PASS")
        cli::cli_h1(
          paste0(
            ppass(),
            " Peaks QC {crayon::bold(x$peaks_analysis$matching_strategy)}: {crayon::green(paste0(prop, '%'))}, {crayon::green(paste('\u03bb =', pur_sc))}. Purity correction: {.value {pur_ch}}."
          )
        )
      else
        cli::cli_h1(
          paste0(
            pfail(),
            " Peaks QC {crayon::bold(x$peaks_analysis$matching_strategy)}: {crayon::red(paste0(prop, '%'))}, {crayon::red(paste('\u03bb =', pur_sc))}. Purity correction: {.value {pur_ch}}."
          )
        )

      cat('\n')

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
        p = format(p * 100, digits = 1)
        p = sprintf("%3s", p)

        cli::cli_alert_info(paste0(
          crayon::blue(xxx$karyotype[1]),
          " ~ n = {n} ({p}%) {clisymbols::symbol$arrow_right} ",
          qc,
          ""
        ))

      }
    }

    # General karyotypes
    if('general' %in% (x$peaks_analysis %>% names))
    {
      gen =  x$peaks_analysis$general$summary %>%
        dplyr::ungroup() %>%
        dplyr::mutate(prop = round(n/sum(n) * 100, 0))

      n_matched =  gen$matched %>% sum
      n_mismatched =  gen$mismatched %>% sum

      cli::cli_h1(
        paste(
          "General peak QC ({.field {sum(gen$n)}} mutations):", ppass(), n_matched, pfail(), n_mismatched, "- epsilon = {.value {x$peaks_analysis$general$params$epsilon}}."
        )
      )
      cat('\n')

      for(i in 1:nrow(gen))
      {
        qc = paste(
          sprintf("%-7s", paste(ppass(), gen$matched[i])),
          sprintf("%-7s", paste(pfail(), gen$mismatched[i]))
        )

        n = sprintf("%-5s", gen$n[i])
        p = sprintf("%3s", gen$prop[i])

        cli::cli_alert_info(
          paste0(
            crayon::blue(gen$karyotype[i]),
            " ~ n = {n} ({p}%) {clisymbols::symbol$arrow_right} ",
            qc,
            ""
          )
        )
      }
    }

    # Subclonal
    if('subclonal' %in% names(x$peaks_analysis))
    {
      nsegs =  x$peaks_analysis$subclonal$expected_peaks$segment_id %>% unique %>% length()

      general_table = x$peaks_analysis$subclonal$summary %>%
        # filter(prop > 0) %>%
        dplyr::group_by(segment_id) %>%
        dplyr::summarise(BR = sum(model == 'branching' & prop > 0), LI = sum(model == 'linear'& prop > 0)) %>%
        dplyr::arrange(dplyr::desc(BR), dplyr::desc(LI))

      certainly_linear = general_table %>% dplyr::filter(BR == 0, LI > 0) %>% nrow()
      certainly_branching = general_table %>% dplyr::filter(BR > 0, LI == 0) %>% nrow()
      ambiguous = general_table %>% dplyr::filter(BR > 0, LI > 0) %>% nrow()

      all_bad = nsegs - (certainly_linear + certainly_branching + ambiguous)

      numcol = function(x)
      {
        if(x == 0) return(crayon::red("{0}"))
        return(paste("{.field {", eval(parse(text=x)), "}}"))
      }

      cli::cli_h1(
        paste(
          "Subclonal peaks QC ({.field {nsegs}} segments, initial state {.field {x$most_prevalent_karyotype}}):",
          crayon::underline("linear"), certainly_linear %>% numcol(),
          crayon::underline("branching"), certainly_branching%>% numcol(),
          crayon::underline("either"), ambiguous%>% numcol(),
          crayon::underline("no support"), all_bad%>% numcol(),
          "- epsilon = {.value {x$peaks_analysis$subclonal$params$epsilon}}."
        )
      )

      puncertain = function()
        "{crayon::bgYellow(crayon::white(\" UNKNOWN \"))}"

      my_print = function(s, b){
       # cat(s)

        hb = x$peaks_analysis$subclonal$summary %>%
          dplyr::filter(segment_id == s) %>%
          dplyr::arrange(dplyr::desc(prop)) %>%
          dplyr::mutate(prop = prop * 100) %>%
          dplyr::mutate(prop = case_when(
            prop == 0 ~ crayon::red('0'),
            prop == 100 ~ crayon::green('100'),
            TRUE ~ crayon::yellow(prop)
          )) %>%
          dplyr::mutate(label = paste0(model_id, " [", prop, "]")) %>%
          dplyr::pull(label) %>%
          paste(collapse = '; ')


        n = x$peaks_analysis$subclonal$summary %>% filter(segment_id == s) %>% filter(row_number() == 1) %>% dplyr::pull(size)
        cl = x$peaks_analysis$subclonal$summary %>% filter(segment_id == s) %>% filter(row_number() == 1) %>% dplyr::pull(clones)
        cl = strsplit(cl, ' ')[[1]]
        cl = paste0(cl[1], ' (', cl[2] %>% as.numeric *100, ') + ', cl[3], ' (', cl[4] %>% as.numeric * 100, ')')
        cl = sprintf("%17s", cl)

        cli::cli_alert_info(
          paste0(
            crayon::blue(sprintf(paste0("%", b, 's'), s)),
            " ~ {n[1]} {crayon::yellow(cl)} : ",
            hb,
            ""
          )
        )
      }

      ml = general_table$segment_id %>% nchar() %>% max

      if(certainly_linear > 0)
      {
        cli::cli_h3(paste(ppass(), "Linear models"))

        general_table %>%
          dplyr::filter(BR == 0, LI > 0) %>%
          dplyr::pull(segment_id) %>%
          lapply(my_print, b = ml)
      }

      if(certainly_branching > 0)
      {
        cli::cli_h3(paste(ppass(), "Branching models"))

        general_table %>%
          dplyr::filter(BR > 0, LI == 0) %>%
          dplyr::pull(segment_id) %>%
          lapply(my_print, b = ml)
      }

      if(ambiguous > 0)
      {
        cli::cli_h3(paste(puncertain(), "Either branching or linear models"))

        general_table %>%
          dplyr::filter(BR > 0, LI > 0) %>%
          dplyr::pull(segment_id) %>%
          lapply(my_print, b = ml)
      }

      if(all_bad > 0)
      {
        cli::cli_h3(paste(pfail(), "Bad models"))

        general_table %>%
          dplyr::filter(BR == 0, LI == 0) %>%
          dplyr::pull(segment_id) %>%
          lapply(my_print, b = ml)
      }

      # nlin = x$peaks_analysis$subclonal$expected_peaks %>%
      #   filter(model == 'linear', matched) %>% nrow()
      # nlin_f = x$peaks_analysis$subclonal$expected_peaks %>%
      #   filter(model == 'linear', !matched) %>% nrow()
      #
      # nbr = x$peaks_analysis$subclonal$expected_peaks %>%
      #   filter(model == 'branching', matched) %>% nrow()
      # nbr_f = x$peaks_analysis$subclonal$expected_peaks %>%
      #   filter(model == 'branching', !matched) %>% nrow()


      # cli::cli_h3(
      #   paste(
      #     "Subclonal peak QC ({.field {nsegs}} segments):",
      #     crayon::bold("linear"),
      #     ppass(), nlin, pfail(), nlin_f,
      #     "~",
      #     crayon::bold("branching"),
      #     ppass(), nbr, pfail(), nbr_f,
      #     "- epsilon = {.value {x$peaks_analysis$subclonal$params$epsilon}}."
      #   )
      # )

    #   S_table = x$peaks_analysis$subclonal$summary
    #   S_table$size = strsplit(S_table$segment_id, split = '\\n') %>%
    #     sapply(function(x) x[[2]])
    #   S_table$clones = strsplit(S_table$segment_id, split = '\\n') %>%
    #     sapply(function(x) x[[3]])
    #   S_table$segment_id = strsplit(S_table$segment_id, split = '\\n') %>%
    #     sapply(function(x) x[[1]]) %>%
    #     strsplit(split = ' ') %>%
    #     sapply(function(x) { paste(x[1], x[2], sep = '@') })
    #
    #   for (s in S_table$segment_id %>% unique())
    #   {
    #     hb = S_table %>%
    #       filter(segment_id == s) %>%
    #       pull(prop) %>%
    #       length() == 1
    #
    #     bp = S_table %>%
    #       filter(segment_id == s) %>%
    #       mutate(prop = prop * 100) %>%
    #       pull(prop) %>%
    #       round(0) %>%
    #       paste0('%')
    #
    #     bm = S_table %>%
    #       filter(segment_id == s) %>%
    #       pull(model)
    #
    #     qc = ifelse(hb,
    #                 crayon::green(paste(bm, bp, collapse = ', ')),
    #                 paste(bm, bp, collapse = ', '))
    #
    #     n = S_table %>% filter(segment_id == s) %>% pull(size)
    #     n = sprintf("%20s", n)
    #     cl = S_table %>% filter(segment_id == s) %>% pull(clones)
    #     cl = sprintf("%17s", cl)
    #
    #     cli::cli_alert_info(
    #       paste0(
    #         crayon::blue(sprintf("%17s", s)),
    #         " ~ {n[1]} {crayon::yellow(cl[1])} {clisymbols::symbol$arrow_right} ",
    #         qc,
    #         ""
    #       )
    #     )
    #
    # }
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

               drv = x$mutations %>% dplyr::filter(karyotype == b)
               if (x %>% has_driver_data())
               {
                 drv = drv %>% dplyr::filter(is_driver) %>% dplyr::pull(driver_label)

                 if (length(drv) > 0)
                   cat('  {', crayon::yellow(paste(drv, collapse = ', ')), '}')
               }

               cat("\n")
             })
}



bar_print_console_scl = function(x, top = nrow(x$cna_subclonal))
{
  max_bars = 10

  L_table =  x$cna_subclonal %>%
    dplyr::mutate(Mb = round(length/1e6, 2), CCF = round(CCF, 2), bars = n/max(n) * max_bars)  %>%
    dplyr::arrange(dplyr::desc(n)) %>%
    dplyr::filter(dplyr::row_number() <= top)

  label_width = 30

  max_n_digits = 1 + (nchar(L_table$n) %>% max)
  max_L_digits = 1 + (nchar(L_table$Mb) %>% max)

  for(i in 1:nrow(L_table))
  {
    cat(sprintf('%*s ', 15, paste(L_table$chr[i],  L_table$from[i], sep = '@')))
    cat(sprintf(' [n = %*s, L = %*s Mb] ',
                max_n_digits,
                L_table$n[i],
                max_L_digits,
                L_table$Mb[i]))

    # str = paste0(L_table$karyotype[i], " (", L_table$CCF[i] ,") - ", L_table$karyotype_2[i],  " (", 1 - L_table$CCF[i] ,") ")
    cat(sprintf('%*s ', 12, paste0(L_table$karyotype[i], " (", L_table$CCF[i] ,")")))
    cat(sprintf('%*s ', 12, paste0(L_table$karyotype_2[i], " (", 1-L_table$CCF[i] ,")")))

    cat(paste(rep("\u25A0", L_table$bars[i]), collapse = ''))
    cat('\n')
  }

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
#' \dontrun{
#' data('example_PCAWG', package = 'CNAqc')
#'
#' plot(example_PCAWG)
#' }
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
