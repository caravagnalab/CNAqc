# #
# # y = readRDS("~/Downloads/dataqc_A4164_1004_PDO.rds")
# # z = readRDS("~/Downloads/dataqc_A3770_1018_PDO.rds")
# # w = readRDS("~/Downloads/dataqc_A3994_117_PDO.rds")
# #
# # x = list(y, z, w)
# #
# # plot_multisample_CNA(x)
# #
# x = readRDS("~/Downloads/all_drivers.rds")
#
# long_format_cna = function(x){
#
#   # LOH
#   LOH = x %>%
#     rowwise() %>%
#     filter(grepl(pattern = '1:0', karyotype)) %>%
#     mutate(type = 'Deletion')
#
#   # AMP
#   AMP = x %>%
#     rowwise() %>% filter((karyotype %in% c("2:0","2:1","2:2") & cluster == "C2")) %>%
#     mutate(type = 'AMP')
#
#   single_copy = x %>%
#     rowwise() %>% filter((karyotype %in% c("2:0","2:1","2:2","1:1") & cluster == "C1")) %>%
#     mutate(type = 'Single_copy')
#
#
#   return(bind_rows(LOH,AMP,single_copy))
#
# }
#
#
# long_format_muts = function(x){
#
#
#   indel = x %>%
#     rowwise() %>%
#     filter(type == "indel")
#
#   snv = x %>%
#     rowwise() %>%
#     filter(type == "SNV")
#   return(bind_rows(indel,snv))
#
# }
#
#
# driver_cna = long_format_cna(x) %>% select(sample,driver_label,type)
#
# driver_type = long_format_muts(x) %>% select(sample,driver_label,type)
#
#
# # prepare oncoprint input
#
# wide_format = bind_rows(driver_cna, driver_type) %>%
#   distinct() %>%
#   group_by(sample, driver_label) %>%
#   mutate(type = paste(type, collapse = ';')) %>%
#   distinct() %>%
#   pivot_wider(names_from = 'driver_label', values_from = "type") %>%
#   ungroup() %>%
#   replace(is.na(.), "") %>%
#   t() %>%
#   data.frame()
#
# colnames(wide_format) = wide_format[1, ]
# wide_format = wide_format[-1, ]
#
#
# ggsci::pal_jama()()
#
# alter_fun_list = list(
#   background = function(x, y, w, h) {
#     grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = grid::gpar(fill = "gainsboro", col = NA))
#   },
#   Deletion = function(x, y, w, h) {
#     grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = grid::gpar(fill = "steelblue", col = NA))
#   },
#   AMP = function(x, y, w, h) {
#     grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = grid::gpar(fill = "indianred3", col = NA))
#   },
#   Single_copy = function(x, y, w, h){
#     grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = grid::gpar(fill = "yellow", col = NA))
#   },
#   SNV = function(x, y, w, h) {
#     grid::grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = grid::gpar(fill = "orange", col = NA))
#   },
#   indel = function(x, y, w, h) {
#     grid::grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = grid::gpar(fill = "forestgreen", col = NA))
#   }
# )
#
# col = c("SNV" = "orange", "indel" = "green", "AMP" = "red", "Deletion" = "blue", "Single_copy" = "yellow")
#
#
#
#
# ComplexHeatmap::oncoPrint(
#   wide_format,
#   get_type = function(x) strsplit(x, ";")[[1]],
#   alter_fun = alter_fun_list,
#   alter_fun_is_vectorized = FALSE,
#   col = col,
#   column_title = "OncoPrint",
#   show_column_names = TRUE,
#   heatmap_legend_param = list(
#     title = "Alternations",
#     at = c("SNV","indel","AMP", "Deletion", "Single_copy"),
#     labels = c("SNV","indel","AMP", "Deletion", "Single_copy")
#   )
# )
