#' get pathways function
#' @param x kegg identifiers
#' @importFrom httr GET message_for_status
get_pathways <- function(x){
  
  name <- x
  
  .strip <- function(str)
  {
    gsub("^\\s+|\\s+$", "", str)
  }
  
  kegg_get <- function(x){
    
    x <- paste(x, collapse="+")
    url <- sprintf("%s/get/%s", getOption("KEGG_REST_URL", "https://rest.kegg.jp"), x)
    response <- GET(url)
    status_report <- tryCatch(message_for_status(GET(url)),
                              http_404 = function(c) "That url doesn't exist",
                              http_403 = function(c) "Authentication required",
                              http_400 = function(c) "Incorrect input from KEGG column",
                              http_500 = function(c) "KEGG server is unavailable"
    )
    content <- .strip(content(response, "text"))
    if (nchar(content) == 0)
      return(status_report)
    return(content)
  }
  
  output <- kegg_get(x)
  
  content <- lapply(output, function(x) strsplit(.strip(x), "\n", fixed=TRUE)[[1]])
  #replace delimeter elements with END_OF_ENTRY to separate entries
  content <- lapply(content, function(x) gsub(x, pattern = "///", replacement = "END_OF_ENTRY"))
  
  #pull out NAME and ORGANISM lines
  content_vect <- content[[1]]
  NAME_content <- content_vect[grep("NAME", x = content_vect)]
  
  NAME_content <- gsub(pattern = "NAME", replacement = "", x = NAME_content)
  
  ORGANISM_content <- content_vect[grep("ORGANISM", content_vect)]
  
  ORGANISM_content <- gsub(pattern = "ORGANISM", replacement = "", x = ORGANISM_content)
  
  #convert to a string
  
  content <- lapply(content, function(x) paste(x, sep = "", collapse = ""))
  #split into character matrix by End of Entry
  
  content <- lapply(content, function(x) str_split(x, "END_OF_ENTRY", simplify = TRUE))
  #remove elements that don't contain REACTION (broken record but needs control flow for each class)
  
  content <- lapply(content, function(x) t(x[,str_detect(x, pattern = "PATHWAY")==TRUE]))
  #convert each column into a vector within a list
  
  #change element names to compound, this will need control flow in the future for each class of cpd, rxn, KO because the word after ENTRY will
  #be different!
  change_names <- function(x){
    
    colnames(x) <- gsub('^.*ENTRY\\s*|\\s*Compound.*$', '', x)
    
    return(x)
    
  }
  content <- lapply(content, change_names)
  content <- lapply(content, function(x) as.list(as.data.frame(x)))
  
  #remove everything but REACTION identifiers (again this will need control flow for each class)
  content <- lapply(content, function(x) lapply(x, function(x) gsub('^.*PATHWAY\\s*|\\s*PATHWAY.*$|MODULE.*$|DISEASE.*$|ENZYME.*$|BRITE.*$|DBLINKS.*$|ATOM.*$|BOND.*$', '', x)))
  
  content <- lapply(content, function(x) lapply(x, function(x) strsplit(x, split = "            ")))
  
  name <- names(content[[1]])
  
  name <- regmatches(name, regexpr(pattern = "EC:.+?(?=])", text = name, perl = TRUE))
  
  name <- paste0(name, "", NAME_content, "", ORGANISM_content)
  
  names(content[[1]]) <- name
  
  pathways <- content[[1]]
  
  pathways[[1]] <- split(pathways[[1]][[1]], c(1:length(pathways[[1]][[1]])))
  
  return(pathways)
}
#' get pathway image function
#' @param path_list list of pathways
#' @param path selected pathway from path_list
#' @importFrom httr GET message_for_status
#' @importFrom magick image_read image_annotate image_ggplot
get_path_image <- function(path_list, path){
  
  pathway <- unlist(path_list[[1]][grep(path, path_list[[1]])])
  
  pathway <- gsub(pattern = " .*", replacement = "", x = path)
  
  kegg_get_image <- function(x){
    
    x <- paste(x, collapse="+")
    url <- sprintf("%s/get/%s/image", getOption("KEGG_REST_URL", "https://rest.kegg.jp"), x)
    response <- image_read(GET(url)$content)
    status_report <- tryCatch(message_for_status(GET(url)),
                              http_404 = function(c) "That url doesn't exist",
                              http_403 = function(c) "Authentication required",
                              http_400 = function(c) "Incorrect input from KEGG column",
                              http_500 = function(c) "KEGG server is unavailable"
    )
    return(response)
    
  }
  
  path_image <- kegg_get_image(pathway)
  
  
  path_image <- image_annotate(path_image, names(path_list), size = 20)
  
  path_image <- image_ggplot(path_image, interpolate = TRUE)
  
  return(path_image)
  
  
  
  
}