

add_context = function(x, reference = BSgenome.Hsapiens.NCBI.GRCh38){
  
  stopifnot(inherits(x, "cnaqc"))
  
  sample = x$sample
  
  snvs = x$mutations 
  
  data <- snvs %>%
    mutate( to = from,
            sample = sample
    ) %>% dplyr::select(sample, chr,from,to,ref,alt,id)
  
  # check that reference is a BSgenome object
  if(is.null(reference)|class(reference)!="BSgenome") {
    stop("The reference genome provided as input needs to be a BSgenome object.")
  }
  
  data$chr = gsub(data$chr,pattern = "chr",replacement = "")
  
  # preprocessing input data
  data <- as.data.frame(data) 
  
  # consider only single nucleotide variants involving (A,C,G,T) bases
  data <- data[which(data[,"from"]==data[,"to"]),,drop=FALSE]
  data <- data[which(as.matrix(data[,"ref"])%in%c("A","C","G","T")),,drop=FALSE]
  data <- data[which(as.matrix(data[,"alt"])%in%c("A","C","G","T")),,drop=FALSE]
  data <- data[,c("sample","chr","from","ref","alt","id"),drop=FALSE]
  colnames(data) <- c("sample","chr","pos","ref","alt","id")
  data <- unique(data)
  data <- data[order(data[,"sample"],data[,"chr"],data[,"pos"]),,drop=FALSE]
  
  # convert data to GRanges
  data <- GRanges(as.character(data$chr),IRanges(start=(as.numeric(data$pos)-1),width=3),
                  ref=DNAStringSet(as.character(data$ref)),alt=DNAStringSet(as.character(data$alt)),sample=as.character(data$sample),id = as.character(data$id))
  
  # check that all chromosomes match reference
  if(length(setdiff(seqnames(data),GenomeInfoDb::seqnames(reference)))>0) {
    warning("Check chromosome names, not all match reference genome.")
  }
  
  # find context for each mutation
  data$context <- getSeq(reference,data)
  
  # check for any mismatch with BSgenome context
  if(any(subseq(data$context,2,2)!=data$ref)) {
    warning("Check reference bases, not all match context.")
  }
  
  # get complements and reverse complements
  data$cref <- complement(data$ref)
  data$calt <- complement(data$alt)
  data$rccontext <- reverseComplement(data$context)
  
  # identify trinucleotides motif
  data$cat <- ifelse(data$ref%in%c("C","T"),paste0(subseq(data$context,1,1),"[",data$ref,">",data$alt,"]",subseq(data$context,3,3)), 
                     paste0(subseq(data$rccontext,1,1),"[",data$cref,">",data$calt,"]",subseq(data$rccontext,3,3)))
  
  sbs = data %>% as_tibble() %>% dplyr::select(id,context) 
  
  x$sbs = full_join(x$mutations %>% filter(type == "SNV"),sbs, by = "id")
  
  x$GRanges = data
  
  return(x)
  
}


get_context_counts = function(x){
  
  if(!"GRanges" %in% names(x)) {
    stop("Run add_context before")
  }
  
  data = x$GRanges
  
  # create 96 trinucleotides mutation categories
  categories_context <- NULL
  categories_alt <- rep(c(rep("C>A",4),rep("C>G",4),rep("C>T",4),rep("T>A",4),rep("T>C",4),rep("T>G",4)),4)
  categories_cat <- NULL
  cont <- 0
  for(i in c("A","C","G","T")) {
    for(j in 1:6) {
      for(k in c("A","C","G","T")) {
        cont <- cont + 1
        categories_context <- c(categories_context,paste0(k,":",i))
        categories_cat <- c(categories_cat,paste0(k,"[",categories_alt[cont],"]",i))
      }
    }
  }
  mutation_categories <- data.table(context=categories_context,alt=categories_alt,cat=categories_cat)
  
  # count number of mutations per sample for each category
  data <- merge(mutation_categories[,.(cat)],data.table(sample=data$sample,cat=data$cat)[,.N,by=.(sample,cat)],by="cat",all=TRUE)
  data <- dcast(data,sample~cat,value.var="N")
  data <- data[!is.na(sample),drop=FALSE]
  data[is.na(data)] <- 0
  
  # make trinucleotides counts matrix
  samples_names <- data$sample
  data <- as.matrix(data[,2:ncol(data),drop=FALSE])
  rownames(data) <- samples_names
  data <- data[sort(rownames(data)),,drop=FALSE]
  data <- data[,sort(colnames(data)),drop=FALSE]
  trinucleotides_counts <- array(0,c(nrow(data),96))
  rownames(trinucleotides_counts) <- rownames(data)
  colnames(trinucleotides_counts) <- sort(as.character(mutation_categories$cat))
  rows_contexts <- rownames(data)
  cols_contexts <- colnames(trinucleotides_counts)[which(colnames(trinucleotides_counts)%in%colnames(data))]
  trinucleotides_counts[rows_contexts,cols_contexts] <- data[rows_contexts,cols_contexts]
  
  # return trinucleotides counts matrix
  return(trinucleotides_counts)
  
}


