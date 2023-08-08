#' UI for uploading and transforming data
#' @param id module id
#' @importFrom shiny NS actionButton HTML numericInput selectInput selectizeInput sliderInput radioButtons downloadButton tagList splitLayout 
#' @importFrom shinyFeedback useShinyFeedback
#' @importFrom shinyWidgets numericInputIcon
upload_and_normalization_UI <- function(id) {
  shinyFeedback::useShinyFeedback()
  ns <- NS(id)
  tagList(
    actionButton(inputId = ns("help_norm"), label = "Help", icon = icon("question-circle", lib = "font-awesome"),style="color: #fff; background-color: #0694bf; border-color: #013747"),
    HTML("<br><br><h4>Upload Data and Create Omu Object</h4><br><br>"),
    splitLayout(fileInput(inputId = ns("upload_metabo"), "Metabolomics Data", accept = c(".csv")),
                    fileInput(inputId = ns("upload_meta"), "Sample Data", accept = c(".csv"))),
    actionButton(inputId = ns("make_list"), "Create Omu Object",style="color: #fff; background-color: #0694bf; border-color: #013747"),
    HTML("<br><br><h4>View Data Distributions</h4><br><br>"),
    selectizeInput(inputId = ns('histo_metabolite'), "Metabolite", choices = NULL, multiple = TRUE),
    splitLayout(selectizeInput(inputId = ns('histo_group'), 'Group By', choices = NULL, multiple = FALSE),
                numericInput(inputId = ns("bins"), label = "Bins", value = 15)),
    radioButtons(inputId = ns("histogram_extension"), "Save Histogram As:",
                 choices = c("png", "pdf", "svg", "pptx"), inline = TRUE),
    downloadButton(ns("download_histogram"), "Download Histogram",style="color: #fff; background-color: #0694bf; border-color: #013747"),
    HTML("<br><br><h4>Normalize Data</h4><br><br>"),
    radioButtons(ns("transform_method"), "Pick Transformation Method:", c("Center by Mean" = "mean_center", "Pareto Scale" = "pareto_scale",
                                                                          "Natural Log" = "ln_t", "glog" = "glog", "Square Root" = "sqrt")),
    actionButton(ns("transform"), "Transform Data",style="color: #fff; background-color: #0694bf; border-color: #013747"),
    HTML("<br><br><h4>Filter Data</h4><br><br>"),
    selectizeInput(inputId = ns("filter_samp"), "Select Samples", choices = NULL, multiple = TRUE),
    actionButton(inputId = ns("remove_samp"), "Remove Samples",style="color: #fff; background-color: #0694bf; border-color: #013747"),
    numericInputIcon(inputId = ns("prev_threshold"), "Select Metabolite Prevalence Across Samples", value = 20, min = 0, max = 100, 
                     icon = list(NULL,icon("percent"))),
    selectizeInput(inputId = ns("filter_metabo"), "Metabolites Below Prevalence Threshold", choices = NULL, multiple = TRUE),
    actionButton(inputId = ns("remove_metabo"), "Remove Metabolites",style="color: #fff; background-color: #0694bf; border-color: #013747")
  )
  
  
}
#'Server for the upload and transform data module
#' @param id module id
#' @importFrom shiny moduleServer observeEvent reactiveValues updateSelectizeInput req updateSelectInput reactive renderUI downloadHandler showModal modalDialog renderPlot
#' @importFrom thematic thematic_shiny
#' @importFrom spsComps shinyCatch
#' @importFrom omu assign_hierarchy
#' @importFrom DT renderDataTable datatable

upload_and_normalization_server <- function(id) {
  
  moduleServer(id, function(input, output, session){
    
    upload_metabo <- reactive({
      inFile <- input$upload_metabo
      
      df <- read_metabo(inFile$datapath)
      #df <- df[,-1]
       return(df)
    })
    
    upload_meta <- reactive({
      inFile <- input$upload_meta
      
      df <- read.csv(inFile$datapath)
      df <- as.data.frame(sapply(df, function(x) gsub(pattern = " ", replacement = ".", x = x)))
      #df <- df[,-1]
      return(df)
    })
    
    #create a reactive data object that forms the base of other server side operations, and include conditions to ensure data is properly formatted
    
    omu_list <- reactiveValues(data = NULL, stats = NULL, rf_list = NULL, kg_data = NULL)
    
    observeEvent(input$make_list,{
      
      #showModal(modalDialog("Creating Omu Object and Assigning Metabolite Metadata", footer=NULL))
      

      metabo <- upload_metabo()
      meta <- upload_meta()
      
      if(any(colnames(metabo)=="KEGG")==FALSE){
        
        shinyCatch(stop("Metabolomics Data missing KEGG column. Fix data and reupload"), blocking_level = "error")
        
      }
      
      if(identical(sort(as.character(colnames(metabo)[unlist(lapply(metabo, is.numeric))])), sort(as.character(meta$Sample))==FALSE)){
        
        shinyCatch(stop("Sample names in metabolomics data and sample data do not match. Fix data and reupload"), blocking_level = "error")
        
      }
      
      if(any(colnames(meta)=="Sample")==FALSE){
        
        shinyCatch(stop("Sample data are missing 'Sample' column. Fix data and reupload"), blocking_level = "error")
        
      }
      
      if(any(colnames(metabo)=="Metabolite")==FALSE){
        
        shinyCatch(stop("Metabolomics data are missing 'Metabolite' column. Fix data and reupload."), blocking_level = "error")
        
      }
      metabo <- assign_hierarchy(count_data = metabo, keep_unknowns = TRUE, identifier = "KEGG")
      
      columns <- c("Class", "Subclass_1", "Subclass_2","Subclass_3", "Subclass_4")
      #eliminate the problem of spaces in names for dynamic color UI
        
      metabo[,names(metabo) %in% columns] <- apply(metabo[,names(metabo) %in% columns], MARGIN = 2, 
                                                   function(x) gsub(pattern = " ", replacement = "_", x = x))
      metabo[,names(metabo) %in% columns] <- apply(metabo[,names(metabo) %in% columns], MARGIN = 2, 
                                                   function(x) gsub(pattern = ",", replacement = "_", x = x))
      
      metabo$Metabolite <- gsub(pattern = "-", replacement = "_", x = metabo$Metabolite)
      
      metabo$Metabolite <- gsub(pattern = "\\(", replacement = "", x = metabo$Metabolite)
      
      metabo$Metabolite <- gsub(pattern = "\\)", replacement = "", x = metabo$Metabolite)
      
      metabo$Metabolite <- gsub(pattern = ",", replacement = "_", x = metabo$Metabolite)
      
      metabo$Metabolite <- gsub(pattern = " ", replacement = "_", x = metabo$Metabolite)
      
      
      #create base_metabo df for volcano plot log2FoldChange
      base_metabo <- metabo
      
      data_list <- list(metabo = metabo, meta = meta,base_metabo = base_metabo)
      omu_list$data <- data_list
      
      
      metabo_prevalence <- metabolite_prevalence(omu_list$data$metabo, threshold = input$prev_threshold)
      
      
      histogram_choices <- omu_list$data$metabo$Metabolite
      updateSelectInput(inputId = "histo_metabolite", choices = histogram_choices, 
                        selected = sample(histogram_choices[grep(pattern = "[a-z]",x = histogram_choices)], 6))
      updateSelectInput(inputId = "histo_group", choices = c(colnames(omu_list$data$meta), "none"), selected = "none")
      updateSelectInput(inputId = "filter_samp", choices = omu_list$data$meta$Sample, selected = "none")
      updateSelectInput(inputId = "filter_metabo", choices = histogram_choices, selected = metabo_prevalence)
      #removeModal()
      
      
    })
    
    #help button pop up text box
    observeEvent(input$help_norm,{
      showModal(modalDialog(
        title = "Help",
        HTML("This is the omu-Shiny module for preparing data for statistical analysis and performing an anova with a tukeys test..<br>
       <br>
       This module has two sections: section 1 is for centering, scaling, and transforming data while section 2 is for
       performing an anova. Histograms are plotted to asses if data are appropriate for the anova, and data can be viewed in the
       data table tab.<br>
       <br>
       When uploading metabolomics data, follow the data format guidelines in the omu vignette, with 2 columns for Metabolite and 
       KEGG compound identifiers, and the rest of the columns as your samples. The Sample Data needs to have a column named 'Sample'
       where the sample names in the sample column match the column headers in the metabolomics data.<br>
       <br>
       Data centering, scaling, and transformation:<br>
       1. Select metabolite(s) to plot and select a sample variable to group histogram bins by.<br>
       <br>
       2. Asses the distribution(s) shown to determine what steps to take to normalize the data. For this its important to understand
       what condition you received the data in. Are these data pre-scaled? Do they need centering? Depending on what has already been done,
       you may only need to log transform the data, or even leave them as is.<br> 
       <br>
       3. The center, pareto, log, and square root action buttons will apply those transformations to the data in real time. To reset
       any of these changes simply click the 'Create Omu Object' button again.<br>" 
        )))
    })
    
    #code for centering, scaling, normalizing data and plotting histogram
    observeEvent(input$transform,{
      
      metabo <- omu_list$data$metabo
      
      transformed_metabo <- switch(input$transform_method,
                                   mean_center = transform_metabolites(metabo, function(x) x - mean(x)),
                                   pareto_scale = transform_metabolites(metabo, function(x) x/sqrt(sd(x))),
                                   ln_t = transform_samples(metabo, log),
                                   glog = transform_samples(metabo,function(x) log(x+1)),
                                   sqrt = transform_samples(metabo, sqrt))
      
      omu_list$data$metabo <- transformed_metabo
      
    })
    
    histogram_plot <- reactive({
      
      #req(is.null(omu_list$data$metabo==FALSE))
      
      
      histogram(metabo = omu_list$data$metabo, meta = omu_list$data$meta, metabolite = input$histo_metabolite, group = input$histo_group, bins = input$bins)
      
    })
    
    
    output$histogram <- renderPlot({
      
      #req(is.null(omu_list$data$metabo==FALSE))
      
      
      #histogram(metabo = omu_list$data$metabo, meta = omu_list$data$meta, metabolite = input$histo_metabolite, group = input$histo_group)
      
      histogram_plot()
    }, bg="transparent", execOnResize = TRUE)
    
    density_p <- reactive({
      
      
      density_plot(metabo = omu_list$data$metabo, meta = omu_list$data$meta, group = input$histo_group)
      
    })
    
    output$density <- renderPlot({
      
      density_p()
      
    })
    
    output$kurtosis_and_skew_table <- DT::renderDataTable({
      
      d <- kurtosis_and_skewness(metabo = omu_list$data$metabo, meta = omu_list$data$meta, metabolite = input$histo_metabolite, group = input$histo_group)
      DT::datatable(data = d, options = list(scrollX = TRUE))
      
    })
    
    observeEvent(input$remove_samp, {
      
      meta <- omu_list$data$meta
      metabo <- omu_list$data$metabo
      base_metabo <- omu_list$data$base_metabo
      
      meta <- meta[!meta$Sample %in% input$filter_samp,]
      metabo <- metabo[,!names(metabo) %in% input$filter_samp]
      base_metabo <- base_metabo[,!names(base_metabo) %in% input$filter_samp]
      
      omu_list$data$metabo <- metabo
      omu_list$data$meta <- meta
      omu_list$data$base_metabo <- base_metabo
    })
    
    observeEvent(input$remove_metabo,{
      
      
      metabo <- omu_list$data$metabo
      metabo <- metabo[!metabo$Metabolite %in% input$filter_metabo,]
      
      base_metabo <- omu_list$data$metabo
      base_metabo <- base_metabo[!base_metabo$Metabolite %in% input$filter_metabo,]
      
      omu_list$data$metabo <- metabo
      omu_list$data$base_metabo <- base_metabo
    })
    
    
    output$download_histogram <- download_plot(input$histogram_extension, histogram_plot())
    
    #a reactive to update the data UI based on which data type  is selected as an input
    data_tables <- reactive({
      
      d <- reactiveValuesToList(omu_list)
      d <- d$data
      
      return(d)
      
    })
    
    output$data_table <- DT::renderDataTable({
      
      d <- data_tables()
      DT::datatable(data = d[["metabo"]], options = list(scrollX = TRUE))
      
    })
    
    output$sample_data <- DT::renderDataTable({
      
      d <- data_tables()
      DT::datatable(data = d[["meta"]], options = list(scrollX = TRUE))
    })
    
    return(omu_list)
  })
  
  
}
#'Produces histogram 
#' @param id module id
#' @importFrom shiny NS plotOutput

histogramOutput <- function(id){
  
  plotOutput(NS(id,"histogram"))
}
#'Produces density plot
#' @param id module id
#' @importFrom shiny NS plotOutput
densityOutput <- function(id){
  
  plotOutput(NS(id, "density"))
  
}
#'Produces kurtosis and skewness table
#' @param id module id
#' @importFrom shiny NS 
#' @importFrom DT dataTableOutput
k_s_table_output <- function(id){
  
  dataTableOutput(NS(id, "kurtosis_and_skew_table"))
  
}
#' metabolite data for data viewer tab
#' @param id module id
#' @importFrom shiny NS 
#' @importFrom DT dataTableOutput
data_viewer_output <- function(id){
  
  dataTableOutput(NS(id, "data_table"))
  
}
#' sample data for data viewer tab
#' @param id module id
#' @importFrom shiny NS 
#' @importFrom DT dataTableOutput
sample_viewer_output <- function(id){
  
  dataTableOutput(NS(id, "sample_data"))
  
}