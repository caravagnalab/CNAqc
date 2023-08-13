#' Maps SNPs bins onto segments 
#'
#' @param snps A dataframe of single nucleotide polymorphism mean BAF and DR over fixed genomic windows (bins) with the following fields:
#'
#' * `chr` chromosome name, e.g., \code{"chr3"}, \code{"chr8"}, \code{"chrX"}, ...
#' * `from` where the bin start, an integer number
#' * `to` where the bin ends, an integer number
#' * `BAF` average B-allele frequency of the bin, a real number in [0,1]
#' * `N.BAF` number if SNPs falling into the bin
#' * `DR` average depth ration of the bin, a real positive number 
#'
#' @param cna A dataframe of allele-specific copy number with the following fields:
#'
#' * `chr` chromosome name, e.g., \code{"chr3"}, \code{"chr8"}, \code{"chrX"}, ...
#' * `from` where the segment start, an integer number
#' * `to` where the segment ends, an integer number
#' * `Major` for the number of copies of the major allele (or A-allele), an integer number
#' * `minor` for the number of copies of the major allele (or B-allele), an integer number
#' * `CCF` an optional cancer cell fraction (CCF) column distinguishing clonal and subclonal segments, a real number in [0,1]
#' * `Major_2` optional for the number of copies of the major allele (or A-allele) in the second clone if present, an integer number
#' * `minor_2` optional for the number of copies of the major allele (or B-allele) in the second clone  if present, an integer number
#'
#' @param subclonal_cna same as cna, but for cna segments with ccf< 1
#'
#' @return the snps object with an extra column, `segment_id` reporting the id of the segment they map on
#' @export
#'
#' @examples
process_snps = function (snps,cna, subclonal_cna = NULL) 
{
  stopifnot(is_tibble(snps) | is.data.frame(snps))
  stopifnot(is_tibble(cna) | is.data.frame(cna))
  #stopifnot(is_tibble(subclonal_cna) | is.data.frame(subclonal_cna))
  
  if (!is.null(subclonal_cna)){
  chromosomes = c(cna$chr, subclonal_cna$chr) %>% unique()
  }else{
    chromosomes = c(cna$chr) %>% unique()
  }
  mapped_snps= lapply(chromosomes, function(c){
    
    snp_c = snps %>% filter(chr==c)
    cna_c = cna %>% filter(chr==c)
    if (!is.null(subclonal_cna)){
      sub_cna_c = subclonal_cna %>% filter(chr==c)
      segments_ids = c(cna_c$segment_id, sub_cna_c$segment_id)
    }else{
      segments_ids = cna_c$segment_id
    }
    new_snp_c = data.frame()
    
    for (i in 1:length(segments_ids)){
      #print('segment : ', i)
      if (segments_ids[i] %in% cna_c$segment_id){
        segment_df = cna_c %>% filter(segment_id == segments_ids[i])
      }else{
        if (!is.null(subclonal_cna) & length(sub_cna_c$segment_id) >0){
          segment_df = sub_cna_c %>% filter(segment_id == segments_ids[i])
        }
      }
      #print(segment_df)
      snp_seg = snp_c %>% filter(from >= segment_df$from, to <= segment_df$to)
      #print(snp_seg)
      new_snp_c = rbind(new_snp_c, snp_seg %>% mutate(segment_id = segments_ids[i]))
      #new_snp_c
    }
    return(new_snp_c)
  })
  
  mapped_snps = Reduce(rbind, mapped_snps)
  
  
  return(mapped_snps)
}

#' Get SNPs dataframe from Sequenza output
#'
#' @param x the result of the following operation:
#' x = load('[sample_name]_sequenza_extract.RData')
#' x = get(x)
#' where the _sequenza_extract.RData file can be found among the outputs of Sequenza
#'
#' @return a list with average BAF and DR either for each bin ($binned) or for each segment ($segment)
#' @export
#'
#' @examples
extract_sequenza_baf_dr = function(x)
{
  r= x
  
  
  segments = r$segments %>%
    Reduce(f = dplyr::bind_rows) %>%
    dplyr::rename(
      chr = chromosome,
      from = start.pos,
      to = end.pos,
      BAF = Bf,
      DR = depth.ratio
    ) %>%
    dplyr::as_tibble()
  
  BAF = r$BAF %>% seq()
  BAF = lapply(BAF, function(k)
    r$BAF[[k]] %>% dplyr::mutate(chr = names(r$BAF)[k]))
  BAF = BAF %>%
    Reduce(f = dplyr::bind_rows) %>%
    dplyr::rename(from = start,
                  to = end,
                  BAF = mean) %>%
    dplyr::select(chr, from, to, BAF, dplyr::everything()) %>%
    dplyr::as_tibble()
  
  DR = r$ratio %>% seq()
  DR = lapply(DR, function(k)
    r$ratio[[k]] %>% dplyr::mutate(chr = names(r$ratio)[k]))
  DR = DR %>%
    Reduce(f = dplyr::bind_rows) %>%
    dplyr::rename(from = start,
                  to = end,
                  DR = mean) %>%
    dplyr::select(chr, from, to, DR, dplyr::everything()) %>%
    dplyr::as_tibble()
  
  return(list(
    segmented = segments,
    binned = dplyr::full_join(
      BAF,
      DR,
      by = c("chr", 'from', 'to'),
      suffix = c(".BAF", ".DR")
    )
  ))
}


#' Returns the segment type of a segment based on CCF and karyotype
#'
#' @param x a CNAqc object
#' @param seg_id the id of the segment 
#'
#' @return the segment type according to CCF and karyotype: 'simple clonal', 'simple subclonal','complex clonal','complex subclonal'
#' @export
#'
#' @examples
get_segment_type = function(x, seg_id){
  seg_type = x$segment_type %>% filter(segment_id == seg_id) %>% pull(segment_type)
  seg_type
}

#' Plot BAF and Depth Ratio of binned SNPs
#'
#' @param x cnaqc object
#' @param what either BAF of DR
#'
#' @return the plot of BAF or DR of binned SNPs along the genome
#' @export
#'
#' @examples
plot_snps = function(x, what='BAF', s= 1)
{
  #BAF_DR = extract_sequenza_baf_dr(x)
  BAF_DR = list(binned = x$snps, segmented = x$cna )
  
  #chrs = BAF_DR$segmented$chr %>% unique()
  chrs = x$cna$chr %>% unique()
  
  BAF_DR$binned= CNAqc:::relative_to_absolute_coordinates(list(reference_genome = 'GRCh38'),
                                                           BAF_DR$binned)
  
  BAF_DR$segmented = CNAqc:::relative_to_absolute_coordinates(list(reference_genome = 'GRCh38'),
                                                              BAF_DR$segmented)
  
  no_x_axis =
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank()
    )
  
  BAF_plot = CNAqc:::blank_genome(ref = 'GRCh38', chromosomes = chrs) +
    ggplot2::geom_hline(
      yintercept = 0.5,
      color = 'indianred3',
      linetype = 'dashed',
      size = .4
    ) +
    ggplot2::geom_segment(
      data = BAF_DR$binned,
      ggplot2::aes(
        x = from,
        xend = to,
        y = BAF,
        yend = BAF
      ),
      size = s,
      colour = 'black'
    ) +
    # ggplot2::geom_segment(
    #   data = BAF_DR$segmented,
    #   ggplot2::aes(
    #     x = from,
    #     xend = to,
    #     y = Bf,
    #     yend = Bf
    #   ),
    #   size = 1,
    #   colour = 'black'
    # ) +
    ggplot2::labs(y = 'BAF')
  # no_x_axis
  
  DR_plot = CNAqc:::blank_genome(ref = 'GRCh38', chromosomes = chrs) +
    ggplot2::geom_hline(
      yintercept = 1,
      color = 'indianred3',
      linetype = 'dashed',
      size = .4
    ) +
    ggplot2::geom_segment(
      data = BAF_DR$binned, #%>% filter(DR <= 2.1),
      ggplot2::aes(
        x = from,
        xend = to,
        y = DR,
        yend = DR
      ),
      size = s,
      colour = 'black'
    ) +
    # ggplot2::geom_segment(
    #   data = BAF_DR$segmented %>% filter(depth.ratio <= 2.1),
    #   ggplot2::aes(
    #     x = from,
    #     xend = to,
    #     y = depth.ratio,
    #     yend = depth.ratio
    #   ),
    #   size = 1,
    #   colour = 'black'
    # ) +
    ggplot2::labs(y = 'Depth ratio')
  # no_x_axis
  
  # cowplot::plot_grid(
  #   BAF_plot ,
  #   DR_plot,
  #   ncol = 1
  # )
  if (what == 'BAF'){
    return(BAF_plot)
  }
  if (what == 'DR'){
    return(DR_plot)
  }
  
  #return(list(BAF_plot, DR_plot))
}

get_segments = function(x, which= c('simple clonal', 'complex clonal', 'simple subclonal', 'complex subclonal')){
  segments = x$segment_type %>% filter(segment_type %in% which) %>% pull(segment_id) #%>% unique()
  segments
}


