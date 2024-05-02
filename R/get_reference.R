get_reference = function(ref)
{

  if(ref %in% c("hg19", "GRCh37"))
    return(CNAqc::chr_coordinates_hg19)

  if(ref %in% c("hg38", "GRCh38"))
    return(CNAqc::chr_coordinates_GRCh38)
  
  if(ref %in% c("mm10", "GRCm38"))
    return(CNAqc::chr_coordinates_GRCm38)
  
  stop("Available references: hg19 (or GRCh37), hg38 (or GRCh38) and mm10 (GRCm38)")
}
