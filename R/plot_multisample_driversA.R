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
# plot_multisample_drivers = function(x)
# {
#   ok_input = sapply(x, function(x)
#     inherits(x, what = 'cnaqc')) %>% all()
#
#   if (!ok_input)
#     stop("Input x must be a list of CNAqc objects!")
#
#   L = x
#   Ln = names(L)
#   if (is.null(Ln)) {
#     cli::cli_alert_warning("The input list is un-named, using default naming scheme Sample*")
#     Ln = paste0("Sample ", 1:length(L))
#
#     names(L) = Ln
#   }
#
#   # Extract all drivers annotations
#   L_flat = lapply(L, Mutations, cna = 'clonal')
#   L_flat = lapply(L_flat %>% names(), function(x) L_flat[[x]] %>% dplyr::mutate(sample = x))
#
#   L_flat = L_flat %>% Reduce(f = bind_rows)
#   if (!("is_driver" %in% colnames(L_flat))) {
#     cli::cli_alert_danger("The input datasets do not contain driver mutations annotated.")
#     return(ggplot())
#   }
#
#   L_flat = L_flat %>% dplyr::filter(is_driver)
#
#   # Pre-process and transform driver annotations
#   x = x %>% dplyr::filter(is_driver)
#
#
#   # Colours adopted
#   CNA_colors = CNAqc:::get_karyotypes_colors(x$karyotype)
#   CNA_colors = c(CNA_colors, `Other` = 'gray')
#
#   mutation_colors = ggsci::pal_lancet()(2)
#   names(mutation_colors) = c("SNV", 'indel')
#
#   colors = c(CNA_colors, mutation_colors)
#
#   # Driver annotations
#   driver_cna = x %>%
#     dplyr::mutate(type = ifelse(karyotype %in% names(CNA_colors), karyotype, "Other")) %>%
#     dplyr::select(sample, driver_label, type)
#
#   driver_muts = x %>%
#     dplyr::select(sample, driver_label, type)
#
#   wide_format = bind_rows(driver_cna, driver_muts) %>%
#     distinct() %>%
#     group_by(sample, driver_label) %>%
#     mutate(type = paste(type, collapse = ';')) %>%
#     distinct() %>%
#     pivot_wider(names_from = 'driver_label', values_from = "type") %>%
#     ungroup() %>%
#     replace(is.na(.), "") %>%
#     t() %>%
#     data.frame()
#
#   colnames(wide_format) = wide_format[1,]
#   wide_format = wide_format[-1,]
#
#   # CNA colors
#   fun_CNA_colors = lapply(CNA_colors,
#                           function(event_color) {
#                             function(x, y, w, h) {
#                               grid::grid.rect(x,
#                                               y,
#                                               w - unit(0.5, "mm"),
#                                               h - unit(0.5, "mm"),
#                                               gp = grid::gpar(fill = event_color, col = NA))
#                             }
#                           })
#
#   # Background colors
#   fun_background_color = list(
#     `background` = function(x, y, w, h) {
#       grid::grid.rect(x,
#                       y,
#                       w - unit(0.5, "mm"),
#                       h - unit(0.5, "mm"),
#                       gp = grid::gpar(fill = "gainsboro", col = NA))
#     }
#   )
#
#   # Mutation colour
#   fun_mutation_colors = lapply(mutation_colors,
#                                function(event_color) {
#                                  function(x, y, w, h) {
#                                    grid::grid.rect(x,
#                                                    y,
#                                                    w - unit(0.5, "mm"),
#                                                    h * 0.33,
#                                                    gp = grid::gpar(fill = event_color, col = NA))
#                                  }
#                                })
#
#   # Alternate list
#   alter_fun_list = append(fun_CNA_colors, fun_background_color) %>% append(fun_mutation_colors)
#
#   # Labels
#   legend_labels = names(alter_fun_list)
#   legend_labels = legend_labels[legend_labels != "background"]
#
#   # pdf(
#   #   paste0("~/Downloads/Driver_statistcs.pdf"),
#   #   width = 7,
#   #   height = 13
#   # )
#   ComplexHeatmap::oncoPrint(
#     wide_format,
#     get_type = function(x)
#       strsplit(x, ";")[[1]],
#     alter_fun = alter_fun_list,
#     alter_fun_is_vectorized = FALSE,
#     col = colors,
#     column_title = "OncoPrint",
#     show_column_names = TRUE,
#     heatmap_legend_param = list(
#       title = "Alternations",
#       at = legend_labels,
#       labels = legend_labels
#     )
#   )
#   dev.off()
# }

plot_multisample_drivers = function(x,
                                    Hugo_Symbol_column = 'gene_symbol',
                                    Variant_Classification_column = 'Variant_Classification')
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

  # maftools mutations object
  maftools_mutations = lapply(1:length(L),
         function(i)
           single_cnaqc_to_maftools(
             x,
             names(x)[i],
             sample_name,
             Hugo_Symbol_column = Hugo_Symbol_column,
             Variant_Classification_column = Variant_Classification_column
           ))



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

  # pdf(
  #   paste0("~/Downloads/Driver_statistcs.pdf"),
  #   width = 7,
  #   height = 13
  # )
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

single_cnaqc_to_maftools = function(x,
                                    sample_name,
                                    Hugo_Symbol_column = 'gene_symbol',
                                    Variant_Classification_column = 'Variant_Classification'
                                    )
{
  # Certain things ar fixed
  mutations = x %>% get_drivers()

  if(nrow(mutations))
  {
    cli::cli_alert_warning("Sample {.field {sample_name}} has no driver mutations ")
  }

  mutations$Chromosome = mutations$chr
  mutations$Start_Position = mutations$from
  mutations$End_Position = mutations$to
  mutations$Reference_Allele = mutations$ref
  mutations$Tumor_Seq_Allele2 = mutations$alt

  # Sample name (from the input named list)
  mutations$Tumor_Sample_Barcode = sample_name

  # Variant_Type -
  # - https://docs.gdc.cancer.gov/Encyclopedia/pages/Variant_Type/
  # - https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/
  # SNP: Single nucleotide polymorphism -- a substitution in one nucleotide
  # DNP: Double nucleotide polymorphism -- a substitution in two consecutive nucleotides
  # TNP: Triple nucleotide polymorphism -- a substitution in three consecutive nucleotides
  # ONP: Oligo-nucleotide polymorphism -- a substitution in more than three consecutive nucleotides
  # INS: Insertion -- the addition of nucleotides
  # DEL: Deletion -- the removal of nucleotides
  Variant_Type = function(mutations, i)
  {
    ref = mutations$ref[i]
    alt = mutations$alt[i]

    nr = nchar(ref)
    na = nchar(alt)

    if(nr == 1 &
       na == 1 &
       ref %in% c("A", "C", "T", "G") &
       alt %in% c("A", "C", "T", "G")
       ) return("SNP")

    if(nr == 2 & na == 2) return("DNP")
    if(nr == 3 & na == 3) return("TNP")
    if(nr >= 4 & na >= 4) return("ONP")

    if(nr > na) return("DEL")
    else return("INS")

    return("UNKNOWN")
  }

  mutations$Variant_Type = sapply(1:nrow(mutations), Variant_Type, mutations = mutations)

  # Genes (from input)
  if(!(Hugo_Symbol_column %in% colnames(mutations)))
    cli::cli_abort("Missing Hugo_Symbol column
                   {.field {Hugo_Symbol_column}} from data")

  if(any(is.na(mutations[[Hugo_Symbol_column]])))
    cli::cli_abort("Some Hugo_Symbol values are NA, remove them.")

  mutations$Hugo_Symbol = mutations[[Hugo_Symbol_column]]

  # Variant_Classification (from input)
  required_Variant_Classification = c(
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "Splice_Site",
    "Translation_Start_Site",
    "Nonsense_Mutation",
    "Nonstop_Mutation",
    "In_Frame_Del",
    "In_Frame_Ins",
    "Missense_Mutation"
  )

  if(!(Variant_Classification_column %in% colnames(mutations)))
    cli::cli_abort("Missing Variant_Classification column
                   {.field {Variant_Classification_column}} from data")

  if(any(is.na(mutations[[Variant_Classification_column]])))
    cli::cli_abort("Some Hugo_Symbol values are NA, remove them.")

  what_misses = mutations[[Variant_Classification_column]]
  what_misses = what_misses[!(what_misses %in% required_Variant_Classification)]

  if(length(what_misses) > 0)
    cli::cli_alert_warning("Some annotated variant classifications
                           are not recognised: {.field {what_misses}}")

  mutations$Variant_Classification = mutations[[Variant_Classification_column]]

  return(mutations)
}

pergene_cnas = function(x, keep_drivers = TRUE, gene_column = 'gene_symbol')
{
  mutations = x %>% Mutations()

  # Genes (from input)
  if(!(gene_column %in% colnames(mutations)))
    cli::cli_abort("Missing gene column
                   {.field {gene_symbol}} from data")

  if(any(is.na(mutations[[gene_symbol]])))
    cli::cli_alert_warning("Some gene_symbol values are NA, will be removed")

  # drivers (if available)
  if(has_driver_data(x) & keep_drivers)
  {
    get_drivers(x)
  }

}
