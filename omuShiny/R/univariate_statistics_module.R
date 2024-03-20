#' UI for univariate statistics
#' @param id module id
#' @importFrom shiny NS actionButton HTML selectizeInput radioButtons downloadButton tagList textInput
#' @importFrom shinyFeedback useShinyFeedback
anova_UI <- function(id) {
  shinyFeedback::useShinyFeedback()
  ns <- NS(id)
  tagList(
    actionButton(inputId = ns("help_anova"), label = "Help", icon = icon("question-circle", lib = "font-awesome"),style="color: #fff; background-color: #0694bf; border-color: #013747"),
    HTML("<br><br><h4>Perform Statistical Tests</h4><br><br>"),
    textInput(inputId = ns("anova_model"), "Enter Model Formula", placeholder = "~ Variable of Interest"),
    radioButtons(ns("stats_test"), "Pick Statistical Model and Post Hoc Test:", c("Anova with Tukeys Post Hoc" = "anova", "Kruskal-Wallis
                                                                                  with Dunn's Post Hoc" = "dunn", "Welch's Anova with Games Powell
                                                                                  Post Hoc" = "welch", "Students T Test" = "students", 
                                                                                  "Welch's T Test" = "welch_t","Wilcox Test" = "wilcox", 
                                                                                  "Paired Students T Test" = "students_p", "Paired Welch's T Test" = "welch_t_p",
                                                                                  "Paired Wilcox Test" = "wilcox_p"
    )),
    actionButton(inputId = ns("run_anova"), "Perform Tests",style="color: #fff; background-color: #0694bf; border-color: #013747"),
    #data viewer
    HTML("<br><br><h4>View Data Tables</h4><br><br>"),
    selectizeInput(inputId = ns("data_viewer"), "View Data", choices = NULL, multiple = FALSE),
    downloadButton(ns("download_anova"), "Download xlsx",style="color: #fff; background-color: #0694bf; border-color: #013747")
  )
}
#'Server for the univariate statistics module
#' @param id module id
#' @importFrom shiny moduleServer observeEvent reactiveValues updateSelectizeInput req updateSelectInput reactive renderUI downloadHandler showModal modalDialog renderPlot
#' @importFrom spsComps shinyCatch
#' @importFrom omu omu_anova omu_summary
#' @importFrom DT renderDataTable datatable
#' @importFrom openxlsx write.xlsx
anova_server <- function(id){
  
  moduleServer(id, function(input, output, session){
    
    
    imputation_check <- function(count_data,metadata,model_characters){
      
      rownames(count_data) <- count_data$Metabolite
      count_data <- count_data[,!names(count_data)%in%c("KEGG")]
      count_data <- split(count_data, f = as.factor(count_data$Metabolite))
      count_data <- lapply(count_data, function(x){x <- x[,!names(x)=="Metabolite"]; return(x)})
      count_data <- lapply(count_data, function(x) as.data.frame(t(x)))
      count_data <- lapply(count_data, function(x){x$Sample <- rownames(x); return(x)})
      count_data <- lapply(count_data, function(x){
        x$term <- metadata[,model_characters][match(x$Sample,metadata$Sample)]; return(x)})
      count_data <- lapply(count_data, function(x){x <- x[-2]; return(x)})
      count_data <- lapply(count_data, function(x){colnames(x)[1] <- "metabolite"; return(x)})
      count_data <- lapply(count_data, function(x) dcast(x, metabolite~term, fun.aggregate = length))
      count_data <- lapply(count_data, function(x){x <- x[-1]; x <- as.data.frame(t(x));
      return(x)})
      count_data <- lapply(count_data, function(x) sapply(x, function(x)x >= 3))
      count_data <- lapply(count_data, function(x) apply(x, MARGIN = 2, all))
      count_data <- lapply(count_data, function(x) any(x)==TRUE)
      count_data <- do.call("rbind", count_data)
      count_data <- any(test_impute[,1])
      return(count_data)
    }
    
    
    anova_model <- reactive({
      req(input$anova_model)
      model <- as.formula(paste("~", input$anova_model))
    
      
    })
    
    test_Factor <- reactive({
      
      Factor <- input$anova_model
      
    })
    
    observeEvent(input$run_anova, {
      if(input$anova_model==""){
        
        shinyCatch(stop("Model term is missing. Type in model term"), blocking_level = "error")
        
      } else{
       model <- anova_model()
        
      }
      req(is.null(omu_list$data)==FALSE)
      #showModal(modalDialog("Running Anova", footer=NULL))
      data <- omu_list$data
      count_data <- data$metabo
      metadata <- data$meta
      #model <- anova_model()
      Factor <- test_Factor()

 
      
      model_characters <- strsplit(gsub("[^[:alnum:] ]", "", as.character(model)[-1]), " +")[[1]]
      
      if(all(model_characters %in% names(metadata))==FALSE){
        
        shinyCatch(stop("One or more model terms do not match column names in metadata. Did you make a typo?"), blocking_level = "error")
        
      }
      
      if(input$stats_test=="dunn"& length(model_characters) > 1){
        
        shinyCatch(stop("Method kruskal can only take a model with one term."), blocking_level = "error")
        
      }
      
      if(input$stats_test=="welch"& length(model_characters) > 1){
        
        shinyCatch(stop("Method welch can only take a model with one term."), blocking_level = "error")
        
      }
      
      if(input$stats_test=="kruskal"&length(levels(as.factor(model_characters[[1]]))) < 3){
        
        shinyCatch(stop("Method kruskal needs at least 3 levels in the model term."), blocking_level = "error")
        
      }
      
      test_impute <- imputation_check(count_data = count_data, metadata = metadata, model_characters = model_characters)
      if(input$stats_test %in% c("welch","students", "students_p", "welch_t", "welch_t_p", "anova") & test_impute==TRUE){
        
        shinyCatch(stop("Repeated values in 3 or more samples for 1 or more metabolites, data might have imputed values. Use wilcox or kruskal test."), blocking_level = "error")
      }
      

     
      omu_list$stats <- switch(input$stats_test, anova = omu_anova(count_data = count_data, metadata = metadata, 
                                                                   model = model, log_transform = FALSE, method = "anova"),
                               dunn = omu_anova(count_data = count_data, metadata = metadata, 
                                                model = model, log_transform = FALSE, method = "kruskal"),
                               welch = omu_anova(count_data = count_data, metadata = metadata, 
                                                 model = model, log_transform = FALSE, method = "welch"),
                               students = list(results = omu_summary(count_data = count_data, metadata = metadata, Factor = Factor,
                                                                     numerator = levels(as.factor(metadata[,Factor]))[1], denominator = levels(as.factor(metadata[,Factor]))[2],
                                                                     test_type = "students")),
                               students_p = list(results =omu_summary(count_data = count_data, metadata = metadata, Factor = Factor,
                                                                      numerator = levels(as.factor(metadata[,Factor]))[1], denominator = levels(as.factor(metadata[,Factor]))[2],
                                                                      test_type = "students", paired = TRUE)),
                               welch_t = list(results =omu_summary(count_data = count_data, metadata = metadata, Factor = Factor,
                                                                   numerator = levels(as.factor(metadata[,Factor]))[1], denominator = levels(as.factor(metadata[,Factor]))[2])),
                               welch_t_p = list(results =omu_summary(count_data = count_data, metadata = metadata, Factor = Factor,
                                                                     numerator = levels(as.factor(metadata[,Factor]))[1], denominator = levels(as.factor(metadata[,Factor]))[2], paired = TRUE)),
                               wilcox = list(results =omu_summary(count_data = count_data, metadata = metadata, Factor = Factor,
                                                                  numerator = levels(as.factor(metadata[,Factor]))[1], denominator = levels(as.factor(metadata[,Factor]))[2], test_type = "mwu")),
                               wilcox_p = list(results =omu_summary(count_data = count_data, metadata = metadata, Factor = Factor,
                                                                    numerator = levels(as.factor(metadata[,Factor]))[1], denominator = levels(as.factor(metadata[,Factor]))[2], test_type = "mwu", paired = TRUE))
      )
      #removeModal()
    })
    
    #a reactive to update the data UI based on which data type  is selected as an input
    data_choice <- reactive({
      
      d <- reactiveValuesToList(omu_list)
      d <- d[["stats"]]
      
      return(d)
      
    })
    #an observation which dynamically updates the data viewer input options when the pick data type input is changed
    observeEvent(data_choice(), {
      
      
      updateSelectInput(inputId = "data_viewer", choices = names(data_choice())) 
      
    })
    
    
    
    output$stats_data_table <- DT::renderDataTable({
      
      d <- data_choice()
      DT::datatable(data = d[[input$data_viewer]], options = list(scrollX = TRUE))
      
    })
    
    output$download_anova <- downloadHandler(
      filename = "omu_data.xlsx",
      
      content = function(file) {
        
        l <- reactiveValuesToList(omu_list)
        l <- c(l$data, l$stats)
        write.xlsx(x = l, file = file)
        
        
      }
    )
    
    observeEvent(input$help_anova,{
      showModal(modalDialog(
        title = "Help",
        HTML("This is the omu-Shiny module for preparing data for statistical analysis and performing an anova with a tukeys test..<br>
       <br>
      Performing an anova:<br> 
               1. Enter a model formula for your anova. For one variable models, simply type the variable name as it is in the Sample Data.
             For multi-variable models, Variable1*Variable2 will give all combinations and comparisons between the two variables. If you want
             the variables modeled separately, do Variable1 + Variable2. For both combinations and separate, Variable1 + Variable2 + Variable1*Variable2.
             This can scale infinitely to however many variables you wish to test. If you have a confounding variable, such as cage effects or paired samples
             in a timecourse, do Variable1 + Confounding_Var, Variable1 + Variable2 + Confounding_Var etc. to subtract that terms variance from the model.<br>
             <br>
               2. To View your data, select data for viewing metabolomics data or sample data, or stats for viewing data from the anova model. Then select the
             specific data table you wish to view.<br> 
               <br> 
               The data can be downloaded to an xlsx spreadsheet where each tab is a table from the omu object. It can also be used in the plotting modules."
        )))
    })
    
    
  })
}

#'Function for viewing stats data
#' @param id module id
#' @importFrom shiny NS 
#' @importFrom DT dataTableOutput
stats_viewer_output <- function(id){
  
  dataTableOutput(NS(id, "stats_data_table"))
  
}