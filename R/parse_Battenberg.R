#' Parse Battenberg calls.
#'
#' @description Returns two tibbles in CNAqc-ready format, one with clonal calls
#' and the other with subclonal calls extracted by parsing the Battenberg format.
#' This function can be used to process raw date released by PCAWG.
#'
#' @param x The input dataframe, required to have the following columns: \code{chr},
#' \code{from}, \code{to}, \code{battenberg_nMaj1_A}, \code{battenberg_nMin1_A},
#' \code{battenberg_nMaj2_A}, \code{battenberg_nMin2_A}, \code{battenberg_frac1_A},
#' and \code{battenberg_frac2_A}.
#'
#' @return Two tibbles in CNAqc-ready format.
#' @export
#'
#' @examples
#' \dontrun{
#' # Load some CSV results from Battenberg
#' x = read.csv(....) %>% parse_Battenberg()
#'
#' # Work with clonal calls (omitting mutations and other parameters here)
#' x = init(cna = x$clonal, ...)
#'
#'
#' # work with subclonal calls (omitting mutations and other parameters here)
#' x = analyze_peaks_subclonal(cna = x$subclonal, ...)
#' }
parse_Battenberg = function(x)
{
  required = c(
    "chr",
    "from",
    "to",
    "battenberg_nMaj1_A",
    "battenberg_nMin1_A",
    "battenberg_nMaj2_A",
    "battenberg_nMin2_A",
    "battenberg_frac1_A",
    "battenberg_frac2_A"
  )

  if (!all(required %in% colnames(x)))
    stop(
      cli::format_error(
        "Missing Battenberg columns in these calls! Required columns are: {.field {required}}"
      )
    )

  all_calls = x %>%
    dplyr::select(chr, from, to, starts_with("Battenberg")) %>%
    dplyr::rename(
      Major = battenberg_nMaj1_A,
      minor = battenberg_nMin1_A,
      Major_2 = battenberg_nMaj2_A,
      minor_2 = battenberg_nMin2_A,
      CCF = battenberg_frac1_A
    ) %>%
    dplyr::select(chr, from, to, Major, minor, Major_2, minor_2, CCF) %>%
    dplyr::filter(!is.na(CCF))

  return(all_calls)

  # clonal_calls = all_calls %>%
  #   dplyr::filter(CCF == 1) %>%
  #   dplyr::rename(Major = Major_1,
  #                 minor = minor_1,
  #                 CCF = CCF_1) %>%
  #   dplyr::select(-Major_2,-minor_2,-CCF_2)
  #
  # subclonal_calls = all_calls %>%
  #   dplyr::filter(!is.na(CCF_2))
  #
  # return(list(clonal = clonal_calls, subclonal = subclonal_calls))
}
