#' get base data function
#' @param base_metabo untransformed metabolite matrix
#' @param meta sample data
#' @param data_stats statistics data
#' @importFrom reshape2 melt
get_base_data <- function(base_metabo, meta,data_stats){
  
  if(any(names(data_stats) %in% "contrast")==TRUE){
    #derive model terms, contrast, numerator, and denominator 
    contrasts <- data_stats$contrast[1]
    contrasts <- unlist(strsplit(contrasts, "-"))
    
    
  }else {
    
    contrasts <- colnames(data_stats)[2:3]
    contrasts <- gsub(".stdev*", "", contrasts)
    
    
  }
  
  numerator <- contrasts[1]
  denominator <- contrasts[2]
  term <- names(meta[,sapply(meta, function(x) any(x %in% contrasts)), drop=FALSE])
  #parse into longform dataframe
  samples_to_remove <- meta[!meta[,term] %in% contrasts,]$Sample
  base_metabo <- base_metabo[,!names(base_metabo) %in% samples_to_remove]
  base_metabo <- reshape2::melt(base_metabo)
  colnames(base_metabo)[9] <- "Abundance"
  colnames(base_metabo)[8] <- "Sample"
  
  base_metabo <- merge(base_metabo, meta, "Sample")
  colnames(base_metabo)[colnames(base_metabo)==term] <- "term"
  return(base_metabo)
  
}




#' caluclate base fold changes
#' @param base_metabo untransformed metabolite matrix
#' @param meta sample data
#' @param data_stats statistics data
#' @importFrom stats aggregate
calc_base_FCs <- function(base_metabo, meta, data_stats){
  
  if(any(names(data_stats) %in% "contrast")==TRUE){
    #derive model terms, contrast, numerator, and denominator 
    contrasts <- data_stats$contrast[1]
    contrasts <- unlist(strsplit(contrasts, "-"))
    
    
  }else {
    
    contrasts <- colnames(data_stats)[2:3]
    contrasts <- gsub(".stdev*", "", contrasts)
    rownames(meta) <- meta$Sample
    index <- sapply(meta, function(x) {all(levels(as.factor(x))%in%contrasts)})
    meta <- meta[,index,drop = FALSE]
    meta$Sample <- rownames(meta)
  }
  
  numerator <- contrasts[1]
  denominator <- contrasts[2]
  
  base_metabo <- base_metabo
  meta <- meta
  data_stats <- data_stats
  
  base_metabo <- get_base_data(base_metabo, meta, data_stats)
  
  base_metabo <- split(base_metabo, f = as.factor(base_metabo$Metabolite))
  #calculate group means
  base_metabo <- lapply(base_metabo, function(x){x <- aggregate(x, Abundance ~ term, mean);
  rownames(x) <- x[,1]
  x <- x[,-1, drop = FALSE]
  x <- as.data.frame(t(x)); 
  return(x)})
  #calculate log2foldchange 
  base_metabo <- lapply(base_metabo, function(x){x$FoldChange <- x[,numerator]/x[,denominator];x$log2FoldChange <- log2(x$FoldChange); return(x)})
  #parse into one data frame and rename columns
  base_metabo <- do.call("rbind", base_metabo)
  colnames(base_metabo)[1] <- paste0("base_mean", "_", colnames(base_metabo)[1])
  colnames(base_metabo)[2] <- paste0("base_mean", "_", colnames(base_metabo)[2])
  colnames(base_metabo)[3] <- paste0("base", "_", colnames(base_metabo)[3])
  colnames(base_metabo)[4] <- paste0("base", "_", colnames(base_metabo)[4])
  
  base_metabo$Metabolite <- rownames(base_metabo)
  
  return(base_metabo)
}

