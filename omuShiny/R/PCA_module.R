#' UI for PCA
#' @param id module id
#' @importFrom shiny NS HTML selectizeInput radioButtons downloadButton tagList splitLayout numericInput checkboxInput 
#' @importFrom shinyFeedback useShinyFeedback
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
                 choices = c("png", "pdf", "svg","eps", "pptx"), inline = TRUE),
    splitLayout(numericInput(ns("figure_height_PCA"), label = "Figure Height(cm)", value = 5),
                numericInput(ns("figure_width_PCA"), label = "Figure Width(cm)", value = 5)),
    downloadButton(ns("download_plot_PCA"), "Save Plot",style="color: #fff; background-color: #0694bf; border-color: #013747"),
    
    
    
  )
  
  
  
}
#'Server for PCA plots
#' @param id module id
#' @importFrom shiny NS moduleServer observeEvent reactiveValuesToList reactiveValues updateSelectizeInput req updateSelectInput reactive renderUI downloadHandler showModal modalDialog renderPlot brushedPoints
#' @importFrom colourpicker colourInput
#' @importFrom ggplot2 ggsave theme theme_bw element_blank element_text element_rect scale_color_manual scale_fill_manual
#' @importFrom omu PCA_plot
#' @importFrom officer read_pptx ph_with add_slide ph_location_type
#' @importFrom DT renderDataTable
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
    
    #dynamic UI that only appears once fill_level inputs are selected
    output$myPanel_PCA <- renderUI({
      
      df<- data()
      pca_var <- input$PCA_variable
      lev <- c(unique(df$meta[,input$PCA_variable])) # sorting so that "things" are unambigious
      cols <- gg_fill_hue(length(lev))
      
      # New IDs "colX1" so that it partly coincide with input$select...
      lapply(seq_along(lev), function(i) {
        colourpicker::colourInput(inputId = NS(id,paste0("col_pca", lev[i])),
                                  label = paste0("Choose Color for ", lev[i]),
                                  value = cbind(cols[i])
        )
        
      })
    })
    
    
    pca <- reactive({
      
      data_plot <- data()
      metabo <- data_plot$metabo
      pca_var <- input$PCA_variable
      cols <- paste0("c(", paste0("input$col_pca", c(unique(data_plot$meta[,input$PCA_variable])), collapse = ", "), ")")
      cols <- eval(parse(text = cols))

      PCA_plot(count_data = data_plot$metabo, metadata = data_plot$meta, variable = input$PCA_variable, color = input$PCA_variable, 
               size = input$PCA_size, label = input$labels, ellipse = input$ellipses) + 
        theme(axis.text = element_text(size = input$PCA_font),axis.title = element_text(size = input$PCA_font),
                                                                                              legend.text = element_text(size = input$PCA_font),
                                                                                              legend.title = element_blank(),
                                                                                              panel.grid = element_blank(),
              panel.border = element_rect(size = input$border_size_PCA), legend.position = "top") +
              scale_color_manual(values = cols) + scale_fill_manual(values = cols)
      
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
                                                                                                 panel.border = element_rect(size = input$border_size_PCA),
                                                                                                 legend.position = "top")
          
          ggsave(file, pca, device = input$extension_PCA,width = input$figure_width_PCA, height = input$figure_height_PCA)
          
        }
        
      }
    )
    
    
  }
  
  
  )
  
  
}
#'Function for pca output
#' @param id module n
#' @importFrom shiny NS plotOutput
pcaOutput <- function(id){
  ns <- NS(id)
  
  plotOutput(ns("PCA_plot"))
  
}
#'Function for pca color output
#' @param id module n
#' @importFrom shiny NS uiOutput
color_output_pca <- function(id){
  
  ns <- NS(id)
  uiOutput(ns("myPanel_PCA"))
  
}