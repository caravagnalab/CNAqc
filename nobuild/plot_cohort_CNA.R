# #
# # y = readRDS("~/Downloads/dataqc_A4164_1004_PDO.rds")
# # z = readRDS("~/Downloads/dataqc_A3770_1018_PDO.rds")
# # w = readRDS("~/Downloads/dataqc_A3994_117_PDO.rds")
# #
# # x = list(y, z, w)
# #
# # plot_multisample_CNA(x)
# #
# x = readRDS("~/Downloads/all_drivers (1).rds")
plot_multisample_drivers = function(x)
{
  ok_input = sapply(x, function(x)
    inherits(x, what = 'cnaqc')) %>% all()

  if (!ok_input)
    stop("Input x must be a list of CNAqc objects!")

  L = x
  Ln = names(L)
  if (is.null(Ln)) {
    cli::cli_alert_warning("The input list is un-named, using default naming scheme Sample*")
    Ln = paste0("Sample ", 1:length(L))

    names(L) = Ln
  }

  # Extract all drivers annotations
  L_flat = lapply(L, Mutations, cna = 'clonal')
  L_flat = lapply(L_flat %>% names(), function(x) L_flat[[x]] %>% dplyr::mutate(sample = x))

  L_flat = L_flat %>% Reduce(f = bind_rows)
  if (!("is_driver" %in% colnames(L_flat))) {
    cli::cli_alert_danger("The input datasets do not contain driver mutations annotated.")
    return(ggplot())
  }

  L_flat = L_flat %>% dplyr::filter(is_driver)

    # Pre-process and transform driver annotations
    x = x %>% dplyr::filter(is_driver)


  # Colours adopted
  CNA_colors = CNAqc:::get_karyotypes_colors(x$karyotype)
  CNA_colors = c(CNA_colors, `Other` = 'gray')

  mutation_colors = ggsci::pal_lancet()(2)
  names(mutation_colors) = c("SNV", 'indel')

  colors = c(CNA_colors, mutation_colors)

  # Driver annotations
  driver_cna = x %>%
    dplyr::mutate(type = ifelse(karyotype %in% names(CNA_colors), karyotype, "Other")) %>%
    dplyr::select(sample, driver_label, type)

  driver_muts = x %>%
    dplyr::select(sample, driver_label, type)

  wide_format = bind_rows(driver_cna, driver_muts) %>%
    distinct() %>%
    group_by(sample, driver_label) %>%
    mutate(type = paste(type, collapse = ';')) %>%
    distinct() %>%
    pivot_wider(names_from = 'driver_label', values_from = "type") %>%
    ungroup() %>%
    replace(is.na(.), "") %>%
    t() %>%
    data.frame()

  colnames(wide_format) = wide_format[1,]
  wide_format = wide_format[-1,]

  # CNA colors
  fun_CNA_colors = lapply(CNA_colors,
                          function(event_color) {
                            function(x, y, w, h) {
                              grid::grid.rect(x,
                                              y,
                                              w - unit(0.5, "mm"),
                                              h - unit(0.5, "mm"),
                                              gp = grid::gpar(fill = event_color, col = NA))
                            }
                          })

  # Background colors
  fun_background_color = list(
    `background` = function(x, y, w, h) {
      grid::grid.rect(x,
                      y,
                      w - unit(0.5, "mm"),
                      h - unit(0.5, "mm"),
                      gp = grid::gpar(fill = "gainsboro", col = NA))
    }
  )

  # Mutation colour
  fun_mutation_colors = lapply(mutation_colors,
                               function(event_color) {
                                 function(x, y, w, h) {
                                   grid::grid.rect(x,
                                                   y,
                                                   w - unit(0.5, "mm"),
                                                   h * 0.33,
                                                   gp = grid::gpar(fill = event_color, col = NA))
                                 }
                               })

  # Alternate list
  alter_fun_list = append(fun_CNA_colors, fun_background_color) %>% append(fun_mutation_colors)

  # Labels
  legend_labels = names(alter_fun_list)
  legend_labels = legend_labels[legend_labels != "background"]


  pdf(
    paste0("~/Downloads/Driver_statistcs.pdf"),
    width = 7,
    height = 13
  )
  ComplexHeatmap::oncoPrint(
    wide_format,
    get_type = function(x)
      strsplit(x, ";")[[1]],
    alter_fun = alter_fun_list,
    alter_fun_is_vectorized = FALSE,
    col = colors,
    column_title = "OncoPrint",
    show_column_names = TRUE,
    heatmap_legend_param = list(
      title = "Alternations",
      at = legend_labels,
      labels = legend_labels
    )
  )
  dev.off()
}
