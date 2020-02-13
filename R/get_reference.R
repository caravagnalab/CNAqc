get_reference = function(ref)
{
  if(ref %in% c("hg19", "GRCh37"))
    return(CNAqc::chr_coordinates_hg19)
  
  if(ref %in% c("GRCh38"))
    return(CNAqc::chr_coordinates_GRCh38)
   
  stop("Available references: hg19 (or GRCh37), and GRCh38")
}