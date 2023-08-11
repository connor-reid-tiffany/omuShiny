#' UI for subset data
#' @param id module id
#' @importFrom shiny NS HTML actionButton selectInput tagList numericInput 
#' @importFrom shinyFeedback useShinyFeedback
subset_data_ui <- function(id){
  
  ns <- NS(id)
  tagList(HTML("<br><br><h4>Select Metabolites</h4><br><br>"),
              actionButton(inputId = ns("start"), "Create Starting Data",style="color: #fff; background-color: #0694bf; border-color: #013747"),
              HTML("<br><br><h5>Subset Metabolites by Metadata</h5><br><br>"),
              selectInput(inputId = ns("metabo_meta"), "Pick Metabolite Metadata Level", 
                          choices = c("Class", "Subclass_1", "Subclass_2", "Subclass_3", "Subclass_4"), multiple = FALSE),
              selectInput(inputId = ns("metabo_meta_level"), "Pick Metabolite Type", choices = NULL, multiple = TRUE),
              actionButton(inputId = ns("subset_meta"), "Subset",style="color: #fff; background-color: #0694bf; border-color: #013747"),
              HTML("<br><br><h5>Subset Metabolites by Random Forest</h5><br><br>"),
              numericInput(inputId = ns("top_n_rf"), label = "Select Top N Important Metabolites", value = 10),
              actionButton(inputId = ns("subset_top_n"), "Subset",style="color: #fff; background-color: #0694bf; border-color: #013747"),
              HTML("<br><br><h5>Subset Metabolites by Statistical Significance</h5><br><br>"),
              selectInput(inputId = ns("stats_d"), label = "Select Statistical Data", choices = NULL, multiple = FALSE),
              numericInput(inputId = ns("sig_level"), label = "Select P Value Threshold", value = 0.05),
              actionButton(inputId = ns("subset_sig"), "Subset by Significance Level",style="color: #fff; background-color: #0694bf; border-color: #013747"))
  
  
}
#' UI for pathway gather
#' @param id module id
#' @importFrom shiny NS HTML actionButton selectInput tagList downloadButton radioButtons
#' @importFrom shinyFeedback useShinyFeedback
pathway_gather_UI <- function(id){
  
  ns <- NS(id)
  tagList(HTML("<br><br><h4>Get Pathway Image</h4><br><br>"),
              HTML("<br><br><h5>Select Gene</h5><br><br>"),
              selectInput(inputId = ns("kg_genes"), "Select Gene", choices = NULL, multiple = FALSE),
              actionButton(inputId = ns("get_pathways"), "Get Pathways ",style="color: #fff; background-color: #0694bf; border-color: #013747"),
              HTML("<br><br><h5>Select and View Pathway</h5><br><br>"),
              selectInput(inputId = ns("select_pathway"), label = "Select Pathway", choices = NULL, multiple = FALSE),
              actionButton(inputId = ns("get_images"), "Get Pathway Image",style="color: #fff; background-color: #0694bf; border-color: #013747"),
              radioButtons(ns("extension_pathway"), "Save As:",
              choices = c("png", "pdf", "svg", "eps"), inline = TRUE),
              downloadButton(ns("download_pathway"), "Save Pathway Image",style="color: #fff; background-color: #0694bf; border-color: #013747"))
  
  
  
}
#' UI for KEGG gather
#' @param id module id
#' @importFrom shiny NS HTML actionButton selectInput tagList downloadButton
#' @importFrom shinyFeedback useShinyFeedback
KEGG_gather_ui <- function(id){
  
  ns <- NS(id)
  tagList(
    HTML("<br><br><h4>Gather Enzyme and Gene Data on Metabolites</h4><br><br>"),
    actionButton(inputId = ns("kegg_gather_p"), "Gather Prokayrotic Data",style="color: #fff; background-color: #0694bf; border-color: #013747"),
    actionButton(inputId = ns("kegg_gather_e"), "Gather Eukaryotic Data",style="color: #fff; background-color: #0694bf; border-color: #013747"),
    HTML("<br><br><h5>Select Enzymes of Interest</h5><br><br>"),
    selectInput(inputId = ns("enzyme_level"), "Select Enzyme Level", choices = c("KO_Class", "KO_Subclass_1", "KO_Subclass_2"), multiple = FALSE),
    selectInput(inputId = ns("enzyme"), "Select Enzymes", choices = NULL, multiple = TRUE),
    actionButton(inputId = ns("subset_enzymes"), "Subset",style="color: #fff; background-color: #0694bf; border-color: #013747"),
    HTML("<br><br><h5>Select Organisms of Interest</h5><br><br>"),
    selectInput(inputId = ns("organism_level"), label = "Select Taxonomic Level", choices = NULL, 
                multiple = FALSE, selected = NULL),
    selectizeInput(inputId = ns("organism"), label = "Select Organisms", choices = NULL, multiple = TRUE),
    actionButton(inputId = ns("subset_genes"), "Subset",style="color: #fff; background-color: #0694bf; border-color: #013747"),
    downloadButton(ns("download_kg"), "Download xlsx",style="color: #fff; background-color: #0694bf; border-color: #013747")
  )
  
  
  
  
}
#'Server for the upload and transform data module
#' @param id module id
#' @importFrom shiny moduleServer observeEvent reactiveValues updateSelectizeInput req updateSelectInput reactive renderUI downloadHandler showModal modalDialog renderPlot
#' @importFrom thematic thematic_shiny
#' @importFrom spsComps shinyCatch
#' @importFrom omu KEGG_gather assign_hierarchy 
#' @importFrom DT renderDataTable datatable
#' @importFrom ggplot2 ggsave
#' @importFrom openxlsx write.xlsx
KEGG_gather_server <- function(id){
  
  moduleServer(id, function(input, output, session){
    
    
    observeEvent(input$start, {
      
      omu_list$kg_data$data <- omu_list$data$metabo
      
      
      
      updateSelectInput(inputId = "metabo_meta_level", choices = unique(omu_list$kg_data$data[,input$metabo_meta]))
      
      
      
      updateSelectInput(inputId = "stats_d", choices = names(omu_list$stats))
      
      
    })
    
    observeEvent(input$subset_meta, {
      
      kg_data <- omu_list$kg_data$data
      
      kg_data <- kg_data[kg_data[,input$metabo_meta] %in% input$metabo_meta_level,]
      
      omu_list$kg_data$data <- kg_data
      
    })
    
    
    observeEvent(input$subset_top_n,{
      
      rf <- omu_list$rf_list
      kg_data <- omu_list$kg_data$data
      importance <- as.data.frame(rf$rf$importance)
      
      if(rf$rf$type=="classification"){
        
        importance <- importance[order(rf$rf$importance[,4],decreasing=TRUE),]
        
      }else if(rf$rf$type=="regression"){
        
        importance <- importance[order(rf$rf$importance[,1],decreasing=TRUE),]
      }
      
      importance$Metabolite <- rownames(importance)
      importance$KEGG <- rf$metabolite_meta$KEGG[match(importance$Metabolite,
                                                       rf$metabolite_meta$Metabolite)]
      
      #importance$OG_names <- omu_list$data$metabo[match(importance$KEGG, omu_list$data$metabo$KEGG)]
      
      importance <- importance[1:input$top_n_rf,]
      
      
      kg_data <- kg_data[kg_data$KEGG %in% importance$KEGG,]
      
      omu_list$kg_data$data <- kg_data
      
    })
    
    observeEvent(input$subset_sig, {
      
      stats <- omu_list$stats
      
      stats_d <- stats[[input$stats_d]]
      
      stats_sig <- stats_d[stats_d$padj <= input$sig_level,]
      
      kg_data <- omu_list$kg_data$data
      
      kg_data <- kg_data[kg_data$KEGG %in% stats_sig$KEGG,]
      
      omu_list$kg_data$data <- kg_data
    })
    
    observeEvent(input$kegg_gather_p,{
      
      
      kg_data <- omu_list$kg_data$data
      
      class(kg_data) <- append(class(kg_data), "cpd")
      
      KG_Enz <- KEGG_gather(kg_data)
      KG_Genes <- KEGG_gather(KG_Enz)
      
      KG_Genes <- assign_hierarchy(KG_Genes, TRUE, identifier = "KO")
      KG_Genes <- assign_hierarchy(KG_Genes, TRUE, identifier = "Prokaryote")
      
      KG_Genes[,sapply(KG_Genes, is.factor)] <- sapply(KG_Genes[,sapply(KG_Genes,is.factor)], as.character)
      
      KG_Genes <- KG_Genes[,!sapply(KG_Genes, is.numeric)]
      #KG_Genes <- KG_Genes[,c(3,6,12,13,14,15,16,17,18,1,2,4,5,7,8,9,10,11)]
      
      omu_list$kg_data$genes <- KG_Genes
      
      
      
    })
    
    observeEvent(input$kegg_gather_e,{
      
      
      kg_data <- omu_list$kg_data$data
      
      class(kg_data) <- append(class(kg_data), "cpd")
      
      KG_Enz <- KEGG_gather(kg_data)
      KG_Genes <- KEGG_gather(KG_Enz)
      
      KG_Genes <- assign_hierarchy(KG_Genes, TRUE, identifier = "KO")
      KG_Genes <- assign_hierarchy(KG_Genes, TRUE, identifier = "Eukaryote")
      
      KG_Genes[,sapply(KG_Genes, is.factor)] <- sapply(KG_Genes[,sapply(KG_Genes,is.factor)], as.character)
      
      KG_Genes <- KG_Genes[,!sapply(KG_Genes, is.numeric)]
      
      #KG_Genes <- KG_Genes[,c(3,6,12,13,14,15,16,17,18,19,1,2,4,5,7,8,9,10,11)]
      
      omu_list$kg_data$genes <- KG_Genes
      
      
      
    })
    6:9
    
    data_genes <- reactive({
      
      data <- reactiveValuesToList(omu_list)
      kg_genes <- data$kg_data$genes
      kg_genes <- kg_genes[,names(kg_genes) %in% c("Kingdom", "Phylum", "Genus.Species", "Common.Name", "Phylum.Class.Family", "Genus", "Species.Strain.Serotype")]
      
      return(kg_genes)
      
    })
    
    
    data_enzymes <- reactive({
      
      data <- reactiveValuesToList(omu_list)
      kg_genes <- data$kg_data$genes
      kg_genes <- kg_genes[,input$enzyme_level]
      
      return(kg_genes)
    })
    
    observeEvent(data_genes(), {
      
      updateSelectInput(inputId = "organism_level", choices = names(data_genes()))
      
      
    })
    
    observeEvent(input$organism_level,{
      req(!is.null(data_genes()))
      updateSelectizeInput(inputId = "organism", choices = unique(data_genes()[,input$organism_level]))
      
    })
    
    observeEvent(data_enzymes(), {
      
      updateSelectInput(inputId = "enzyme", choices = unique(data_enzymes()))
      
    })
    
    observeEvent(input$subset_enzymes,{
      
      KG_Genes <- omu_list$kg_data$genes
      
      KG_Genes <- KG_Genes[KG_Genes[,input$enzyme_level] %in% input$enzyme,]
      
      omu_list$kg_data$genes <- KG_Genes
      
    })
    
    observeEvent(input$subset_genes,{
      
      KG_Genes <- omu_list$kg_data$genes
      
      KG_Genes <- KG_Genes[KG_Genes[,input$organism_level] %in% input$organism,]
      
      omu_list$kg_data$genes <- KG_Genes
      
      data_g <- omu_list$kg_data$genes
      
      genes <- data_g$Genes
      
      updateSelectInput(inputId = "kg_genes", choices = unique(genes))
      
    })
    
    observeEvent(input$get_pathways,{
      
      #subset data to observation with input$kg_genes
      data_g <- omu_list$kg_data$genes
      
      data_g <- as.data.frame(data_g[data_g$Genes==input$kg_genes,])
      
      #paste org and gene values with :
      gene <- paste0(data_g$Org, ":", data_g$Genes)
      
      #remove any gene name by targetting parenthesis and characters after with regex
      gene <- gsub(pattern = "\\(.*", replacement = "", x = gene)
      
      #get pathways
      pathways <- get_pathways(gene)
      
      omu_list$kg_data$pathways <- pathways
    })
    
    data_pathways <- reactive({
      
      data <- reactiveValuesToList(omu_list)
      kg_pathways<- data$kg_data$pathways
      
      
      return(kg_pathways)
    })
    
    observeEvent(data_pathways(), {
      
      data <- data_pathways()
      data_unlist <- as.character(unlist(data[[1]]))
      
      updateSelectInput(inputId = "select_pathway", choices = data_unlist)
      
    })
    
    observeEvent(input$get_images,{
      
      path_image <- get_path_image(path_list = data_pathways(), path = input$select_pathway)
      
      omu_list$kg_data$path_image <- path_image
    })
    
    output$path_image <- renderPlot(res = 96, width = 1200, height = 1200,expr = {
      
      d <- reactiveValuesToList(omu_list)
      i <- d$kg_data$path_image
      
      return(i)
      
    })
    
    #a reactive to update the data UI based on which data type  is selected as an input
    kegg_tables <- reactive({
      
      d <- reactiveValuesToList(omu_list)
      d <- d$kg_data
      
      return(d)
      
    })
    
    output$kg_data_table <- DT::renderDataTable({
      
      d <- kegg_tables()
      DT::datatable(data = d[["data"]], options = list(scrollX = TRUE))
      
    })
    
    output$kg_genes_table <- DT::renderDataTable({
      
      d <- kegg_tables()
      DT::datatable(data = d[["genes"]], options = list(scrollX = TRUE))
    })
    
    output$download_kg <- downloadHandler(
      filename = "kg_data.xlsx",
      
      content = function(file) {
        
        l <- reactiveValuesToList(omu_list)
        l <- c(l$kg_data)
        write.xlsx(x = l, file = file, append = TRUE)
        
        
      }
    )
    
    output$download_pathway <- downloadHandler(
      filename = function() {
        paste("Pathway", input$extension_pathway, sep = ".")
      },
      
      content = function (file) {
    
        ggsave(file, omu_list$kg_data$path_image, device = input$extension_pathway)
          
        }
        
      
    )
    
    
    
  })
  
  
}
#'Function for kg data table
#' @param id module n
#' @importFrom shiny NS 
#' @importFrom DT dataTableOutput
KG_viewer_output <- function(id){
  
  dataTableOutput(NS(id, "kg_data_table"))
  
}
#'Function for genes table
#' @param id module n
#' @importFrom shiny NS 
#' @importFrom DT dataTableOutput
genes_viewer_output <- function(id){
  
  dataTableOutput(NS(id, "kg_genes_table"))
  
}
#'Function for synthetic pathway image
#' @param id module n
#' @importFrom shiny NS plotOutput
pathway_image_output <- function(id){
  
  plotOutput(NS(id, "path_image"))
  
}