advanced_phasing = function(x, cutoff_n = 50)
{
  karyotypes = x$n_karyotype[x$n_karyotype >= cutoff_n] %>% names()

  phasing =  lapply(karyotypes,
                    function(karyotype) {

                      cat("\n")
                      cli::cli_h1("Advanced phasing for {.field {karyotype}}")
                      cat("\n")

                      cli::cli_alert("Signal deconvolution")

                      # x_k = x %>% subset_by_segment_karyotype(karyotype)
                      x_k = x$mutations %>% dplyr::filter(karyotype == !!karyotype)

                      n_ploidy = strsplit(karyotype, ':')[[1]] %>% as.numeric() %>% max
                      n_ploidy = 1 + n_ploidy

                      x_fit = BMix::bmixfit(
                        data = x_k %>% dplyr::select(NV, DP) %>% as.data.frame(),
                        K.Binomials = 1:n_ploidy,
                        samples = 2,
                        silent = TRUE
                      )

                      # cli::cli_h2("Detected {.field {x_fit$K[1]}} peaks for {.field {karyotype}}")

                      # Expectations
                      peaks = expectations_generalised(p = x$purity, karyotype = karyotype)

                      for (i in 1:nrow(peaks))
                        peaks$cluster[i] = abs(x_fit$B.params - peaks$peak[i]) %>% which.min() %>% names()

                      peaks$p_bin = x_fit$B.params[peaks$cluster]
                      peaks$offset = x_fit$B.params[peaks$cluster] - peaks$peak

                      peaks = peaks %>%
                        dplyr::arrange(cluster, abs(offset)) %>%
                        dplyr::group_by(cluster) %>%
                        dplyr::filter(dplyr::row_number() == 1) %>%
                        dplyr::ungroup()

                      p_v = peaks$multiplicity
                      names(p_v) = peaks$cluster
                      
                      # Add potential tail cluster with minimum multiplicity
                      missing_label = setdiff(x_fit$labels %>% unique, names(p_v))
                      if(length(missing_label) > 0)
                      {
                        n_missing = missing_label %>% length()
                        
                        np_v = rep(min(p_v), n_missing)
                        names(np_v) = setdiff(x_fit$labels %>% unique, names(p_v))
                        
                        p_v = c(p_v, np_v)
                      }
                      
                      print(peaks %>% dplyr::select(multiplicity, peak, cluster, offset))

                      # 762 2351
                      # xx_k = x_k %>% select(VEP.SYMBOL, is_driver, VAF) 
                      # xx_k$multiplicity = x_fit$labels
                      # 
                      # xx_k %>%
                      #   dplyr::mutate(
                      #     # multiplicity = x_fit$labels[dplyr::row_number()],
                      #     # multiplicity = p_v[multiplicity],
                      #     # CCF = ccf_adjustment_fun(
                      #     #   v = VAF,
                      #     #   m = strsplit(karyotype, ':')[[1]][2] %>% as.numeric(),
                      #     #   M = strsplit(karyotype, ':')[[1]][1] %>% as.numeric(),
                      #     #   p = x$purity,
                      #     #   mut.allele = multiplicity
                      #     # )
                      #   ) %>% 
                      #   filter(is_driver)
                      
                      x_k$multiplicity = x_fit$labels
                      
                      x_k %>%
                        dplyr::mutate(
                          # multiplicity = x_fit$labels[dplyr::row_number()],
                          multiplicity = p_v[multiplicity],
                          CCF = ccf_adjustment_fun(
                            v = VAF,
                            m = strsplit(karyotype, ':')[[1]][2] %>% as.numeric(),
                            M = strsplit(karyotype, ':')[[1]][1] %>% as.numeric(),
                            p = x$purity,
                            mut.allele = multiplicity
                          )
                          )

                    })

  phasing = Reduce(dplyr::bind_rows, phasing)

  # melted_phasing = phasing %>%
  #   dplyr::select(VAF,CCF, karyotype, driver_label, multiplicity, is_driver) %>%
  #   dplyr::mutate(id = dplyr::row_number()) %>%
  #   reshape2::melt(id = c('id', "karyotype", "driver_label", "multiplicity", 'is_driver'))
  #
  # melted_phasing %>%
  #   ggplot2::ggplot() +
  #   ggplot2::geom_histogram(
  #     ggplot2::aes(value, fill = multiplicity %>% paste),
  #     binwidth = 0.01) +
  #   ggplot2::facet_grid(karyotype~variable, scales = 'free') +
  #   CNAqc:::my_ggplot_theme() +
  #   ggplot2:: guides(fill = ggplot2::guide_legend("Multiplicity")) +
  #   ggplot2::labs(title = "Multiplicity phasing") +
  #   ggsci::scale_fill_jama()

  plot_m = phasing %>%
    ggplot2::ggplot() +
    ggplot2::geom_histogram(
      ggplot2::aes(VAF, fill = multiplicity %>% paste),
      binwidth = 0.01) +
    ggplot2::facet_wrap(~karyotype, scales = 'free') +
    CNAqc:::my_ggplot_theme() +
    ggplot2:: guides(fill = ggplot2::guide_legend("Multiplicity")) +
    ggplot2::labs(title = "Multiplicity phasing", y = "counts") +
    ggsci::scale_fill_jama()

  if(nrow(phasing %>% filter(is_driver)) > 0)
    plot_m = annotate_phased_drivers(x, plot_m, phasing)

  plot_ccf = phasing %>%
    ggplot2::ggplot() +
    ggplot2::geom_histogram(
      ggplot2::aes(CCF, fill = multiplicity %>% paste),
      binwidth = 0.01) +
    ggplot2::facet_wrap(~karyotype, scales = 'free') +
    CNAqc:::my_ggplot_theme() +
    ggplot2::guides(fill = ggplot2::guide_legend("Multiplicity")) +
    ggplot2::labs(title = "CCF per segment",y = "counts") +
    ggsci::scale_fill_jama()

  if(nrow(phasing %>% filter(is_driver)) > 0)
    plot_ccf = annotate_phased_drivers(x, plot_ccf, phasing)

  colors = get_karyotypes_colors(c("1:0", "2:0", '1:1', '2:1', '2:2', '0:0'))
  missing = setdiff(phasing$karyotype %>% unique(), colors %>% names)
  extras = ggsci::pal_jco()(missing %>% length)
  names(extras) = missing
  colors = c(colors, extras)

  plot_ccf_overall = phasing %>%
    ggplot2::ggplot() +
    ggplot2::geom_histogram(
      ggplot2::aes(CCF, fill = karyotype), binwidth = 0.01
      ) +
    CNAqc:::my_ggplot_theme() +
    ggplot2::guides(fill = ggplot2::guide_legend("Segment")) +
    ggplot2::labs(title = "CCF",y = "counts") +
    ggplot2::scale_fill_manual(values = colors)

  if(nrow(phasing %>% filter(is_driver)) > 0)
  {
    # CCF special
    ymax = ggplot_build(plot_ccf_overall)$layout$panel_params[[1]]$y.range[2]

    drivers = phasing %>%
      dplyr::filter(is_driver) %>%
      dplyr::mutate(
          driver_label = ifelse(
            nchar(driver_label) > 10,
            paste0(substr(driver_label, 0, 10), '..'),
            driver_label
          )
        )

    plot_ccf_overall = plot_ccf_overall +
        ggrepel::geom_label_repel(
          data = drivers,
          ggplot2::aes(
            x = CCF,
            y =  1,
            label = driver_label,
            fill = karyotype %>% paste
          ),
          color = 'white',
          ylim = c(ymax * .7, ymax * .9),
          size = 2,
          nudge_x = 0,
          show.legend = FALSE,
          segment.size = 0.25,
          segment.linetype = 2,
          segment.curvature = 1,
          segment.ncp = 1,
          segment.square = TRUE,
          segment.inflect = TRUE,
          segment.color = 'black'
        )
  }


  figure = cowplot::plot_grid(
    plot_m,
    plot_ccf,
    plot_ccf_overall,
    nrow = 1,
    align = 'h',
    axis = 'tb'
  )

  return(
    list(
      phasing = phasing,
      plot = figure
    )
  )

}

annotate_phased_drivers = function(x, plot, phasing)
{
  # VAF/DP/NV etc
  ggb = ggplot_build(plot)
  facet_p = ggb$layout$layout %>% dplyr::mutate(y_max = NA)

  for(i in 1:nrow(facet_p))
    facet_p$y_max[i] = ggb$layout$panel_params[[i]]$y.range[2]

  drivers = phasing %>%
    dplyr::filter(is_driver) %>%
    dplyr::select(VAF, driver_label, karyotype, multiplicity) %>%
    dplyr::left_join(
      facet_p %>% dplyr::select(karyotype, y_max),
      by = c('karyotype')
    ) %>%
    as_tibble() %>% 
    dplyr::group_split(karyotype)

  for(i in 1:length(drivers))
  {
    d = drivers[[i]] %>%
      dplyr::mutate(
        driver_label = ifelse(
          nchar(driver_label) > 10,
          paste0(substr(driver_label, 0, 10), '..'),
          driver_label
        )
      )

    plot = plot +
      ggrepel::geom_label_repel(
        data = d,
        ggplot2::aes(
          x = VAF,
          y =  1,
          label = driver_label,
          fill = multiplicity %>% paste
        ),
        color = 'white',
        ylim = c(d$y_max[1] * .7, d$y_max[1] * .9),
        size = 2,
        nudge_x = 0,
        show.legend = FALSE,
        segment.size = 0.25,
        segment.linetype = 2,
        segment.curvature = 1,
        segment.ncp = 1,
        segment.square = TRUE,
        segment.inflect = TRUE,
        segment.color = 'black'
      )
  }

  return(plot)
}

