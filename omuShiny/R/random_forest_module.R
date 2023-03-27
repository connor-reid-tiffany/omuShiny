random_forest_ui <- function(id){
  
  ns <- NS(id)
  tagList(
    HTML("<br><br><h4>Perform Random Forest</h4><br><br>"),
    textInput(inputId = ns("rf_model"), "Enter Model Formula", placeholder = "Variable of Interest ~."),
    selectInput(inputId = ns("remove_NA"), "Remove Unidentified Metabolites", choices = c("TRUE", "FALSE"), selected = "FALSE"),
    numericInput(inputId = ns("train_prop"), "Choose Percentage of Samples for Training Data", value = 80, min = 1, max = 100),
    actionButton(inputId = ns("rf"), "Perform Random Forest",style="color: #fff; background-color: #0694bf; border-color: #013747"),
    HTML("<br><br><h4>Adjust Plot Parameters</h4><br><br>"),
    selectInput(inputId = ns("rf_meta"),"Select Metabolite Metadata for Color", choices =  c("Class", "Subclass_1", "Subclass_2","Subclass_3", "Subclass_4")),
    numericInput(inputId = ns("rf_n"), "Choose n Metabolites to Plot", value = 10, min = 1, max = 50),
    splitLayout( checkboxInput(inputId = ns("rf_labels"), label = "Toggle Point Labels", value = FALSE),
                 checkboxInput(inputId = ns("rf_ellipses"), label = "Toggle CI Ellipses", value = FALSE)
    ),
    splitLayout(numericInput(inputId = ns("rf_point_size"), label = "Point Size", value = 2.5),
                numericInput(inputId = ns("rf_font"), label = "Font Size", value = 14)
    ),
   radioButtons(ns("extension_rf"), "Save As:",
                 choices = c("png", "pdf", "svg", "pptx"), inline = TRUE),
    splitLayout( numericInput(ns("figure_height_rf"), label = "Figure Height(cm)", value = 5, width = "70%"),
                 numericInput(ns("figure_width_rf"), label = "Figure Width(cm)", value = 5, width = "70%")),
    downloadButton(ns("download_plot_rf"), "Save Plot",style="color: #fff; background-color: #0694bf; border-color: #013747")
    
    
  )
  
}

random_forest_server <- function(id){
  
  moduleServer(id, function(input, output, session){
    
    data <- reactive({
      
      data <- reactiveValuesToList(omu_list)
      data <- omu_list$data
      
      return(data)
      
    })
    
    model <- reactive({
      
      mod <- as.formula(paste(input$rf_model, "~."))
      
    })
    
    proportion <- reactive({
      
      train <- input$train_prop
      
      test <- 100 - train
      
      vect <- c(train, test)
      
      return(vect)
      
    })
    
    observeEvent(input$rf,{
      
      df <- data()
      
      meta <- df$meta
      meta <- as.data.frame(sapply(meta, as.character))
      meta[] <- lapply(meta, as.factor)
      
      if(input$remove_NA==FALSE){
        
        metabo <- df$metabo
        
      }else if(input$remove_NA==TRUE){
        
        metabo <- df$metabo
        metabo <- metabo[metabo$KEGG!="NA",]
        metabo <- metabo[!is.na(metabo$KEGG),]
        
      }
      
      omu_list$rf_list <- omu::random_forest(count_data = metabo, metadata = meta, model = model(), training_proportion = proportion(), n_tree = 500)
      updateSelectInput(inputId = "rf_group", choices = colnames(meta[,names(meta)!="Sample", drop = FALSE]))
      
    })
    
    rf_data <- reactive({
      
      data <- reactiveValuesToList(omu_list)
      rf_data <- data$rf_list
      return(rf_data)
    })
    
    variable_importance <- reactive({
      req(rf_data())
      rf <- rf_data()
      plot <- plot_variable_importance(rf_list = rf, color = input$rf_meta, n_metabolites = input$rf_n)+ 
        scale_x_discrete(label = function(x){ x <- gsub(pattern = "ID_", replacement = "", x = x); x <- stringr::str_trunc(x, 17); return(x)}) + labs(color = NULL, size = "Mean Decrease\nAccuracy") 
      plot <- plot + theme(plot.tag.position = c(0,1),panel.grid = element_blank(), panel.border = element_rect(size = 1.5), 
                           legend.position = "right", legend.box = "vertical",legend.justification = "left", 
                           legend.box.just = "left",axis.text = element_text(size = input$rf_font), axis.title.x = element_text(size = input$rf_font), 
                           legend.text = element_text(size = input$rf_font-5), legend.title = element_text(size = input$rf_font-5), 
                           legend.spacing.x = unit(-0.6, 'mm'), legend.margin = ggplot2::margin(0,3,0,0), legend.box.spacing = unit(0,"pt")
                           ,legend.box.background = element_blank(), legend.key = element_rect(color = NA,fill=NA))+
                    guides(color = guide_legend(override.aes = list(size=1, linetype = 0, fill = NA, legend.position = "top")))
      return(plot)
    })
    
    output$var_imp_plot <- renderPlot(res = 96,{
      
      variable_importance()
      
    }, bg = "transparent")
    
    
    rf_pca <- reactive({
      
      req(rf_data())
      rf <- rf_data()
      req(rf$rf$type=="classification")
      
      plot <- plot_rf_PCA(rf_list = rf, color = input$rf_model, size = input$rf_point_size, label = input$rf_labels, ellipse = input$rf_ellipses)
      
      plot <- plot + theme(panel.grid = element_blank(), panel.border = element_rect(size = 1.5), legend.position = "top",
                                        axis.text = element_text(size = input$rf_font-1), axis.title = element_text(size = input$rf_font), 
                                        legend.text = element_text(size = input$rf_font),legend.title = element_blank()) 
      return(plot)
      
    })
    
    output$rf_PCA <- renderPlot(res = 96,{
      
      rf_pca()
      
    }, bg = "transparent")
    
    output$download_plot_rf <- downloadHandler(
      filename = function() {
        paste("rf_plot", input$extension_rf, sep = ".")
      },
      
      content = function (file) {
        if (grepl(".pptx", file)==TRUE){
          doc <-  read_pptx() 
          doc <- add_slide(doc, layout =  'Title and Content', master = 'Office Theme')
          doc <- ph_with(doc, value = arrangeGrob(variable_importance(), rf_pca(), nrow = 2), location = ph_location_type(type = "body")) 
          print(doc, file)
          
        }else {
          var_imp <- variable_importance() + theme(plot.background = element_rect(color = "white"))
          rf <- rf_pca() + theme(plot.background = element_rect(color = "white"))
          rf_plots <- cowplot::plot_grid(var_imp, rf, nrow = 2,align = "hv",axis = "r")
          
          ggsave(file,rf_plots, device = input$extension_rf,width = input$figure_width_rf, height = input$figure_height_rf)
          
        }
        
      }
    )
    
  }) 
}

varimpOutput <- function(id){
  ns <- NS(id)
  
  plotOutput(ns("var_imp_plot"))
  
}

rfPCAOutput <- function(id){
  ns <- NS(id)
  
  plotOutput(ns("rf_PCA"))
  
}