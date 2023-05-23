has_MAF_annotations = function(x)
{
  grepl("MAF.", colnames(x$mutations)) %>% sum() > 3
}

