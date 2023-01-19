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
    checkboxInput(inputId = ns("rf_labels"), label = "Toggle Point Labels", value = FALSE),
    checkboxInput(inputId = ns("rf_ellipses"), label = "Toggle CI Ellipses", value = FALSE),
    sliderInput(inputId = ns("rf_point_size"), label = "Point Size", min = 0.5, max = 5, value = 2.5),
    sliderInput(inputId = ns("rf_font"), label = "Font Size", min = 5, max = 50, value = 14)
    
    
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
      plot <- plot_variable_importance(rf_list = rf, color = input$rf_meta, n_metabolites = input$rf_n) + ggtitle("Variable Importance") + theme(title = element_text(size = 20))
      plot <- plot + xlab("Metabolite") + theme(panel.grid = element_blank(), panel.border = element_rect(size = 1.5), legend.position = "right", 
                                        axis.text = element_text(size = input$rf_font), axis.title = element_text(size = input$rf_font), 
                                        legend.text = element_text(size = input$rf_font), legend.title = element_text(size = input$rf_font))
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
                                        axis.text = element_text(size = input$rf_font), axis.title = element_text(size = input$rf_font), 
                                        legend.text = element_text(size = input$rf_font),legend.title = element_text(size = input$rf_font))
      return(plot)
      
    })
    
    output$rf_PCA <- renderPlot(res = 96,{
      
      rf_pca()
      
    }, bg = "transparent")
    
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