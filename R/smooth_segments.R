#' Smooth copy number segments.
#'
#' @description
#' This functions join segments that have the same Major and minor alleles
#' (absolute copy number), and that are split by at most a certain number
#' of nucleotides \code{Delta}. The pre-smoothing copy number segments are
#' retained in the output computation.
#'
#' @param x An object of class \code{cnaqc}, created by the \code{init} function.
#' @param maximum_distance The \code{Delta} maximum of nucleotides to allow
#' to join two equivalent copy number segments.
#'
#' @return An object of class \code{cnaqc}, created by the \code{init} function.
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna,example_dataset_CNAqc$purity)
#'
#' x = smooth_segments(x)
#' plot_smoothing(x)
smooth_segments = function(x, maximum_distance = 1e6)
{
  segments = x$cna
  smoothed_segments = NULL

  ncnacl = sum(segments$CCF == 1)
  ncnasbcl = sum(segments$CCF < 1)

  # Subclonal CNA calls -- raise informative warning
  if(ncnasbcl > 0)
    cli::boxx("Subclonal CNAs detected in the dataset, those segments will NOT be removed for the smoothing process.
Remove them before calling 'CNAqc::init' if you want to smoothe only clonal segments.", col = 'red')


  for(chr in unique(segments$chr))
  {
    chr_segments = segments %>% filter(chr == !!chr)

    # Special case, one segment
    if(nrow(chr_segments) == 1) {
      smoothed_segments = bind_rows(smoothed_segments, chr_segments)
      next
    }

    # cat('\n')
    # cli::cli_alert_info("Smoothing {.field {chr}}: {.value {nrow(chr_segments)}} segments.")
    cat("Smoothing", crayon::blue(chr), "with", crayon::red(nrow(chr_segments)), "segments: ")

    # General case: read every segment, start from 1, index tracks where we
    # start merging, another points moves ahead to detect segments ot merge
    index = 1

    repeat{
      # The template is the segment with index "index",
      # from and karyotypes will always be from that
      template = chr_segments[index, ]

      # j moves ahead as far as possible, starting from index
      # when we stop moving j, we merge all from index to j (inclusive)
      j = index

      repeat{

        # Stop if end of list
        if(j == nrow(chr_segments)) break

        separation = (chr_segments$from[j + 1] - chr_segments$to[j]) < maximum_distance
        minor_match = chr_segments$minor[j + 1] == chr_segments$minor[j]
        Major_match = chr_segments$Major[j + 1] == chr_segments$Major[j]

        # We move to the next if we can merge, and stop otherwise
        if(separation & minor_match & Major_match) j = j + 1
        else break
      }

      # cli::cli_alert("Smoothed segments from {.field {index}} to {.value {j}}")
      if(j > index) cat(paste0("[", index, '-',j, '] '))

      # ending here, j is the index of the last segment we merge (inclusive)
      template$to = chr_segments$to[j]
      template$length = template$to - template$from

      smoothed_segments = bind_rows(smoothed_segments, template)

      # Stop if we added the last possible segment
      if(j == nrow(chr_segments)) break

      # Otherwise the new index is j + 1
      index = j + 1
    }

    cat("\n")

  }

  cat('\n')
  cli::cli_alert_success("Smoothed from {.value {nrow(segments)}} to {.value {nrow(smoothed_segments)}} segments with {.value {maximum_distance}} gap (bases).")
  cli::cli_alert_info("Creating a new CNAqc object. The old object will be retained in the $before_smoothing field.")

  clonal_CNA = smoothed_segments %>% dplyr::select(-segment_id, -n)
  subclonal_CNA = x$cna_subclona %>% dplyr::select(-segment_id, -n, -analysed)

  cna = bind_rows(clonal_CNA, subclonal_CNA)%>% dplyr::select(-starts_with('karyotype'), -mutations)

  # Clean up the new segments table,
  x_new = CNAqc::init(x$snvs,
                      cna,
                      purity = x$purity,
                      ref = x$reference_genome)
  x_new$before_smoothing = x

  return(x_new)
}

# # short implementation..
# ...
#
# # Gap shorter than maximum_distance
# deltas = chr_segments$from[2:nrow(chr_segments)] - chr_segments$to[1:(nrow(chr_segments) - 1)]
# deltas = deltas < maximum_distance
#
# # Same copy numer
# minors = sapply(2:nrow(chr_segments), function(x) chr_segments$minor[x] == chr_segments$minor[x-1])
# Majors = sapply(2:nrow(chr_segments), function(x) chr_segments$Major[x] == chr_segments$Major[x-1])
#
# chr_segments$can_merge[2:nrow(chr_segments)] = deltas & minors & Majors

