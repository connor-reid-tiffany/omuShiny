download_plot <- function(extension, plot){
  
  downloadHandler(
    filename = function() {
      paste("plot", extension, sep = ".")
    },
    
    content = function (file) {
      if (grepl(".pptx", file)==TRUE){
        
        doc <-  read_pptx() 
        doc <- add_slide(doc, layout =  'Title and Content', master = 'Office Theme')
        doc <- ph_with(doc, value = plot, location = ph_location_type(type = "body")) 
        print(doc, file)
        
      }else {
        
        ggsave(file, plot, device = extension)
        
      }
      
    }
  )
  
  
  
}