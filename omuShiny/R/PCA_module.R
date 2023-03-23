PCA_ui <- function(id){
  
  ns <- NS(id)
  tagList(
    HTML("<br><br><h4>Principle Component Analysis</h4><br><br>"),
    selectizeInput(inputId = ns("PCA_variable"), "Pick Variable to Group by", choices = NULL, multiple = FALSE),
    splitLayout(checkboxInput(inputId = ns("labels"), label = "Toggle Point Labels", value = FALSE),
                checkboxInput(inputId = ns("ellipses"), label = "Toggle CI Ellipses", value = FALSE)),
    splitLayout(numericInput(inputId = ns("PCA_size"), label = "Point Size", value = 2.5),
                numericInput(inputId = ns("PCA_font"), label = "Font Size", value = 14),
                numericInput(ns("border_size_PCA"), "Border Size", min = 0.5, max = 5, value = 1.5)),
    radioButtons(ns("extension_PCA"), "Save As:",
                 choices = c("png", "pdf", "svg", "pptx"), inline = TRUE),
    splitLayout(numericInput(ns("figure_height_PCA"), label = "Figure Height(cm)", value = 5),
                numericInput(ns("figure_width_PCA"), label = "Figure Width(cm)", value = 5)),
    downloadButton(ns("download_plot_PCA"), "Save Plot",style="color: #fff; background-color: #0694bf; border-color: #013747"),
    
    
    
  )
  
  
  
}

PCA_server <- function(id){
  
  moduleServer(id, function(input, output, session){
    
    data <- reactive({
      
      data <- reactiveValuesToList(omu_list)
      data <- omu_list$data
      
      return(data)
      
    })
    
    observeEvent(data(), {
      
      df <- data()
      meta <- df$meta
      
      updateSelectizeInput(inputId = "PCA_variable", choices = colnames(meta[,names(meta)!="Sample",drop = FALSE]), selected = NULL)
      
    })
    
    
    pca <- reactive({
      
      data_plot <- data()
      metabo <- data_plot$metabo

      PCA_plot(count_data = data_plot$metabo, metadata = data_plot$meta, variable = input$PCA_variable, color = input$PCA_variable, 
               size = input$PCA_size, label = input$labels, ellipse = input$ellipses) + 
        theme(axis.text = element_text(size = input$PCA_font),axis.title = element_text(size = input$PCA_font),
                                                                                              legend.text = element_text(size = input$PCA_font),
                                                                                              legend.title = element_blank(),
                                                                                              panel.grid = element_blank(),
              panel.border = element_rect(size = input$border_size_PCA))
      
    })
    
    output$PCA_plot <- renderPlot(res = 96,{
      
      pca()
      
    }, bg = "transparent")
    
    #download handler for plot, the boolean after content is required for pptx to be an option, as it is not covered by ggsave
    output$download_plot_PCA <- downloadHandler(
      filename = function() {
        paste("PCA_plot", input$extension_PCA, sep = ".")
      },
      
      content = function (file) {
        if (grepl(".pptx", file)==TRUE){
          
          doc <-  read_pptx() 
          doc <- add_slide(doc, layout =  'Title and Content', master = 'Office Theme')
          doc <- ph_with(doc, value = pca(), location = ph_location_type(type = "body")) 
          print(doc, file)
          
        }else {
          
          pca <- pca() + theme_bw() + theme(panel.grid = element_blank(), panel.background = element_rect(fill = "white"), 
                                            plot.background = element_rect(fill = "white")) + 
            theme(axis.text = element_text(size = input$PCA_font),axis.title = element_text(size = input$PCA_font),
                                                                                                 legend.text = element_text(size = input$PCA_font),
                                                                                                 legend.title = element_blank(),
                                                                                                 panel.grid = element_blank(),
                                                                                                 panel.border = element_rect(size = input$border_size_PCA))
          
          ggsave(file, pca, device = input$extension_PCA,width = input$figure_width_PCA, height = input$figure_height_PCA)
          
        }
        
      }
    )
    
    
  }
  
  
  )
  
  
}

pcaOutput <- function(id){
  ns <- NS(id)
  
  plotOutput(ns("PCA_plot"))
  
}