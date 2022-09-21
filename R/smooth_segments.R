#' Smooth simple clonal CNAs.
#'
#' @description
#' This function joins simple clonal CNAs that have the same Major and minor alleles,
#' and that are zeparate by at most a certain number
#' of nucleotides. The pre-smoothing copy number segments are
#' retained in the output computation.
#'
#' @param x A CNAqc object.
#' @param maximum_distance The maximum number of nucleotides to  join two equivalent
#' segments.
#'
#' @return A CNAqc object with smoothed segments, which retains also the pre-smooting
#' dataset.
#' @export
#'
#' @examples
#' data('example_dataset_CNAqc', package = 'CNAqc')
#' x = init(mutations = example_dataset_CNAqc$mutations, cna = example_dataset_CNAqc$cna, purity = example_dataset_CNAqc$purity)
#'
#' # Before smoothing
#' print(x)
#' x = smooth_segments(x)
#'
#' # After smoothing
#' plot_smoothing(x)
smooth_segments = function(x, maximum_distance = 1e6)
{
  segments = x$cna
  smoothed_segments = NULL

  # Subclonal CNA calls -- raise informative warning
  if(x$n_cna_subclonal > 0)
    cli::boxx("Subclonal CNAs will not be smoothed!", col = 'red')

  for(chr in unique(segments$chr))
  {
    chr_segments = segments %>% filter(chr == !!chr)
    chr_segments_subc = x$cna_subclonal %>% filter(chr == !!chr)

    chr_segments = chr_segments %>% bind_rows(chr_segments_subc) %>% arrange(from)

    # Special case, one segment
    if(nrow(chr_segments) == 1) {
      smoothed_segments = bind_rows(smoothed_segments, chr_segments)
      next
    }

    # cat('\n')
    # cli::cli_alert_info("Smoothing {.field {chr}}: {.value {nrow(chr_segments)}} segments.")
    # cat("Smoothing", crayon::blue(chr), "with", crayon::red(nrow(chr_segments)), "segments: ")

    cat("\u2192", crayon::blue(chr),  crayon::red(nrow(chr_segments)), "-")

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

        if("Major_2" %in% colnames(chr_segments))
        {
          is_subclonal_next = !is.na(chr_segments$Major_2[j+1]) & !is.na(chr_segments$minor_2[j+1])
          is_subclonal_this = !is.na(chr_segments$Major_2[j]) & !is.na(chr_segments$minor_2[j])

          if(is_subclonal_next | is_subclonal_this) break
        }

        separation = (chr_segments$from[j + 1] - chr_segments$to[j]) < maximum_distance
        minor_match = chr_segments$minor[j + 1] == chr_segments$minor[j]
        Major_match = chr_segments$Major[j + 1] == chr_segments$Major[j]

        # We move to the next if we can merge, and stop otherwise
        if(separation & minor_match & Major_match) j = j + 1
        else break
      }

      # cli::cli_alert("Smoothed segments from {.field {index}} to {.value {j}}")
      # if(j > index) cat(paste0("[", index, '-',j, '] '))

      # ending here, j is the index of the last segment we merge (inclusive)
      template$to = chr_segments$to[j]
      template$length = template$to - template$from

      smoothed_segments = bind_rows(smoothed_segments, template)

      # Stop if we added the last possible segment
      if(j == nrow(chr_segments)) break

      # Otherwise the new index is j + 1
      index = j + 1
    }

    cat(crayon::green((smoothed_segments$chr %>% table())[chr]), "@")


    cat("\n")
  }

  cat('\n')
  cli::cli_alert_success("Smoothed from {.value {nrow(segments)}} to {.value {nrow(smoothed_segments)}} segments with {.value {maximum_distance}} gap (bases).")
  cli::cli_alert_info("Creating a new CNAqc object. The old object will be retained in the $before_smoothing field.")

  if("mutations" %in% colnames(smoothed_segments)) smoothed_segments = smoothed_segments %>% dplyr::select(-mutations)
  if(grepl('karyotype', colnames(smoothed_segments)) %>% any) smoothed_segments = smoothed_segments %>% dplyr::select(-mutations)
  if("segment_id" %in% colnames(smoothed_segments)) smoothed_segments = smoothed_segments %>% dplyr::select(-segment_id)
  if("n" %in% colnames(smoothed_segments)) smoothed_segments = smoothed_segments %>% dplyr::select(-n)

  # clonal_CNA = smoothed_segments %>% dplyr::select(-segment_id, -n)

  # Extract subclonal CNAs
  # subclonal_CNA = NULL
  # if(!is.null(x$cna_subclonal) & nrow(x$cna_subclonal) > 0)
  # {
  #   subclonal_CNA = x$cna_subclonal %>% dplyr::select(-segment_id, -n, -analysed)
  #   cna = bind_rows(clonal_CNA, subclonal_CNA) %>% dplyr::select(-starts_with('karyotype'), -mutations)
  # }
  # else
  #   cna = clonal_CNA %>% dplyr::select(-starts_with('karyotype'))

  # Extract mutations mapped to subclonal CNAs
  # subclonal_mutations_CNA = NULL
  # if(!is.null(x$cna_subclonal) & nrow(x$cna_subclonal) > 0)
  #   subclonal_mutations_CNA = Reduce(bind_rows, x$cna_subclonal$mutations)
  #
  # muattions = bind_rows(x$mutations, subclonal_mutations_CNA)


  # Clean up the new segments table,
  x_new = init(
    mutations = x %>% Mutations(),
    cna = smoothed_segments,
    purity = x$purity,
    ref = x$reference_genome,
    sample = x$sample)

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

