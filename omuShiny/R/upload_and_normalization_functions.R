#' Histogram function
#' @param metabo metabolite matrix
#' @param meta sample data
#' @param metabolite selected metabolites
#' @param bins number of bins to plot
#' @param group independent variable to group by
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot geom_histogram aes facet_wrap ylab theme element_blank element_text element_rect
histogram <- function(metabo, meta, metabolite, bins = 30, group="none"){
  
  #subset to metabolite of choice
  metabo_sub <- metabo[metabo$Metabolite %in% metabolite,]
  
  #melt
  metabo_sub_melt <- reshape2::melt(metabo_sub, value.name = "Metabolite_Intensity")
  
  #merge melt with metadata
  colnames(metabo_sub_melt)[8] <- "Sample"
  metabo_sub_melt <- merge(metabo_sub_melt, meta, by = "Sample")
  
  
  #plot
  
  if(group=="none"){
    plot <- ggplot(metabo_sub_melt,aes(x = Metabolite_Intensity)) + geom_histogram(color = "black",fill = "white", size = 2, bins = bins) + 
      facet_wrap(.~Metabolite, scales = "free") + 
      ylab("Number of Samples per Bin") +
      xlab("Abundance") +
      theme_bw() +
        theme(panel.grid = element_blank(), axis.text = element_text(size = 16), axis.text.x = element_text(angle = 90), axis.title = element_text(size = 20), strip.text = element_text(size = 10),
                         legend.title = element_blank(), legend.text = element_text(size = 14),         panel.background = element_rect(fill = "transparent",colour = NA),
              plot.background = element_rect(fill = "transparent",colour = NA))
  }else if(group!="none"){

    plot <- ggplot(metabo_sub_melt,aes(x = Metabolite_Intensity)) + geom_histogram(data = metabo_sub_melt,aes_string(color = group), fill = "white", size = 2, bins = bins)+ 
      facet_wrap(.~Metabolite, scales = "free") + 
      ylab("Number of Samples per Bin") +
      xlab("Abundance") +
      theme_bw() +
      theme(panel.grid = element_blank(), axis.text = element_text(size = 16),axis.text.x = element_text(angle = 90), axis.title = element_text(size = 20), strip.text = element_text(size = 10),
                         legend.title = element_blank(), legend.text = element_text(size = 14),        panel.background = element_rect(fill = "transparent",colour = NA),
            plot.background = element_rect(fill = "transparent",colour = NA))
  }
  return(plot)
  
  
}
#' function to determine metabolite prevalence across samples
#' @param metabo metabolite matrix
#' @param threshold user supplied threshold for subsetting metabolites

metabolite_prevalence <- function(metabo, threshold){
  
  rownames(metabo) <- metabo$Metabolite
  metabo <- metabo[,sapply(metabo, is.numeric)]
  
  metabo <- as.data.frame(t(metabo))
  
  prev_df <- data.frame(Metabolite = colnames(metabo),Prevalence = sapply(metabo, function(x) length(x > 0)/ length(x)))
  
  prev_df <- prev_df[which(prev_df$Prevalence < threshold/100),]
  
  return(prev_df)
}
#' density plot function
#' @param metabo metabolite matrix
#' @param meta sample data
#' @param group independent variable to group by
#' @importFrom ggplot2 ggplot geom_density aes facet_wrap ylab theme element_blank element_text element_rect

density_plot <- function(metabo, meta, group){
  
  rownames(metabo) <- metabo$Metabolite
  metabo <- metabo[,sapply(metabo, is.numeric)]
  
  metabo_sum <- as.data.frame(sapply(metabo, sum))
  
  metabo_sum$Sample <- rownames(metabo_sum)
  colnames(metabo_sum)[1] <- "Total_Metabolite_Abundance"
  metabo_sum <- merge(metabo_sum, meta, "Sample")
  
  if(group=="none"){
    
    plot <- ggplot(data = metabo_sum, aes(x = Total_Metabolite_Abundance)) + geom_density(size =2,alpha = 1/3,color = "midnightblue", fill = "dodgerblue2") +
      theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18), legend.position = "none", axis.title.x = element_blank())
    
  }else if(group!="none"){
    
    plot <- ggplot(data = metabo_sum, aes(x = Total_Metabolite_Abundance, color = metabo_sum[,group], fill = metabo_sum[,group])) + 
      geom_density(size = 2,alpha = 1/3) +
      theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18), legend.position = "none", axis.title.x = element_blank() )
    
  }
  
  return(plot)
}
#' kurtosis and skewness function
#' @param metabo metabolite matrix
#' @param meta sample data
#' @param group independent variable to group by
#' @param metabolite user selected metabolites
#' @importFrom reshape2 melt
#' @importFrom stats aggregate p.adjust
#' @importFrom car leveneTest
#' @importFrom e1071 kurtosis skewness
kurtosis_and_skewness <- function(metabo, meta, metabolite, group){
  
  #subset to metabolite of choice
  metabo_sub <- metabo[metabo$Metabolite %in% metabolite,]
  
  #melt
  metabo_sub_melt <- reshape2::melt(metabo_sub, value.name = "Metabolite_Intensity")
  
  #merge melt with metadata
  colnames(metabo_sub_melt)[which(names(metabo_sub_melt) == "variable")] <- "Sample"
  metabo_sub_melt <- merge(metabo_sub_melt, meta, by = "Sample")
  
  
  
  if(group=="none"){
    
    
    kurtosis_tab <- aggregate(Metabolite_Intensity ~ Metabolite, metabo_sub_melt, kurtosis)
    colnames(kurtosis_tab)[2] <- "kurtosis"
    
    skewness_tab <- aggregate(Metabolite_Intensity ~ Metabolite, metabo_sub_melt, skewness)
    colnames(skewness_tab)[2] <- "skewness"
    
    kurt_and_skew <- merge(kurtosis_tab, skewness_tab, by = "Metabolite")
    
  }else if(group!="none"){
    
    metabo_sub_melt$Metabolite_and_Group <- paste0(metabo_sub_melt$Metabolite, "_", metabo_sub_melt[,group])
    
    kurtosis_tab <- aggregate(Metabolite_Intensity ~ Metabolite_and_Group, metabo_sub_melt, kurtosis)
    colnames(kurtosis_tab)[2] <- "kurtosis"
    
    skewness_tab <- aggregate(Metabolite_Intensity ~ Metabolite_and_Group, metabo_sub_melt, skewness)
    colnames(skewness_tab)[2] <- "skewness"
    
    metabo_sub_melt_list <- split(metabo_sub_melt, f = factor(metabo_sub_melt$Metabolite))
    
    test_levenes <- lapply(metabo_sub_melt_list, function(x) leveneTest(Metabolite_Intensity ~ Metabolite_and_Group, x))
    
    test_levenes <- lapply(test_levenes, function(x) x[[3]][1])
    
    test_levenes <- as.data.frame(do.call("rbind", test_levenes))
    
    colnames(test_levenes)[1] <- "Levene_Test_padj"
    
    test_levenes$Levene_Test_padj <- p.adjust(test_levenes$Levene_Test_padj, method = "BH")
    
    test_levenes$Metabolite <- rownames(test_levenes)
    
    kurt_and_skew <- merge(kurtosis_tab, skewness_tab, by = "Metabolite_and_Group")
    
    kurt_and_skew$Metabolite <- gsub(pattern = "_[^_]+$", replacement = "", x = kurt_and_skew$Metabolite_and_Group, perl = TRUE)
    
    kurt_and_skew <- merge(kurt_and_skew, test_levenes, by = "Metabolite")
    
  }
  
  return(kurt_and_skew)
  
}
