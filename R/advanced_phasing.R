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
                      x_k = x$snvs %>% dplyr::filter(karyotype == !!karyotype)

                      n_ploidy = strsplit(karyotype, ':')[[1]] %>% as.numeric() %>% max

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

                      print(peaks %>% dplyr::select(multiplicity, peak, cluster, offset))

                      x_k %>%
                        dplyr::mutate(
                          multiplicity = x_fit$labels[dplyr::row_number()],
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

  plot_m = phasing %>%
    ggplot() +
    geom_histogram(aes(VAF, fill = multiplicity %>% paste), binwidth = 0.01) +
    facet_wrap(~karyotype, scales = 'free') +
    CNAqc:::my_ggplot_theme() +
    guides(fill = guide_legend("Multiplicity")) +
    labs(title = "Multiplicity phasing") +
    ggsci::scale_fill_jama()

  plot_ccf = phasing %>%
    ggplot() +
    geom_histogram(aes(CCF, fill = multiplicity %>% paste), binwidth = 0.01) +
    facet_wrap(~karyotype, scales = 'free') +
    CNAqc:::my_ggplot_theme() +
    guides(fill = guide_legend("Multiplicity")) +
    labs(title = "CCF per segment") +
    ggsci::scale_fill_jama()

  plot_ccf_overall = phasing %>%
    ggplot() +
    geom_histogram(aes(CCF, fill = karyotype), binwidth = 0.01) +
    CNAqc:::my_ggplot_theme() +
    guides(fill = guide_legend("Segment")) +
    labs(title = "CCF") +
    ggsci::scale_fill_lancet()

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

