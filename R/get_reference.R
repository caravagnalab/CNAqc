get_reference = function(ref, data)
{

  if(ref %in% c("hg19", "GRCh37"))
    return(CNAqc::chr_coordinates_hg19)

  if(ref %in% c("hg38", "GRCh38"))
    return(CNAqc::chr_coordinates_GRCh38)
  
  if(ref %in% c("mm10", "GRCm38"))
    return(CNAqc::chr_coordinates_GRCm38)
  
  if(!ref %in% c("hg19", "GRCh37", "hg38", "GRCh38", "mm10", "GRCm38")) {
    check_reference(data = data)
    return(data)
  }
  
  # stop("Available references: hg19 (or GRCh37), hg38 (or GRCh38) and mm10 (GRCm38)")
}


check_reference = function(data) {
  if(any(!colnames(data) %in% c("chr", "length", "from", "to", "centromerStart", "centromerEnd"))) {
    cli::cli_abort("Missing correct colnames in input data")
  }
}
