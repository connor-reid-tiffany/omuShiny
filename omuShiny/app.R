#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyFeedback)
library(shinyWidgets)
library(omu)
library(tidyverse)
library(openxlsx)
library(DT)
library(officer)
library(thematic)
library(colourpicker)
library(e1071)
library(car)
library(httr)
library(magick)
library(spsComps)
library(ggpubr)
library(ggrepel)
library(gridExtra)
library(stringr)

#' The shiny app
#' @param ... placeholder
#' @importFrom shiny fluidPage titlePanel tags radioButtons sidebarLayout sidebarPanel tabsetPanel tabPanel mainPanel observe session shinyApp verticalLayout conditionalPanel
#' @importFrom bslib bs_theme bs_theme_update
myApp <- function(...){
  my_theme <- bslib::bs_theme(bootswatch = "flatly")
ui <- fluidPage(
  titlePanel("omu-Shiny: Analyze and Plot your Metabolomics Data"),
  tags$style(type="text/css",
             ".shiny-output-error { visibility: hidden; }",
             ".shiny-output-error:before { visibility: hidden; }"),
  theme = my_theme,
  radioButtons("current_theme", "App Theme:", c("Light" = "flatly", "Dark" = "darkly"), inline = TRUE),
  sidebarLayout(
    sidebarPanel(tags$head(tags$style(type="text/css", "
             #loadmessage {
               position: fixed;
               top: 0px;
               left: 0px;
               width: 100%;
               padding: 5px 0px 5px 0px;
               text-align: center;
               font-weight: bold;
               font-size: 100%;
               color: #000000;
               background-color: #CCFF66;
               z-index: 105;
             }
          ")),conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                               tags$div("Loading...",id="loadmessage")),
      tabsetPanel(
        tabPanel("Upload and Transform Data", upload_and_normalization_UI("myModule1")),
        tabPanel("Machine Learning", tabsetPanel(type = "pills",
                                                 tabPanel("PCA",PCA_ui("myModule1"),color_output_pca("myModule1")), tabPanel("Random Forest", random_forest_ui("myModule1")))),
        tabPanel("Univariate Statistics", anova_UI("myModule1")),
        tabPanel("Univariate Plots", tabsetPanel(
          tabPanel("Volcano Plot", volcano_ui("myModule1"),color_output("myModule1")),
          tabPanel("Dotplot", dotplot_ui("myModule1"), color_output_dot("myModule1")),
          tabPanel("Boxplot", boxplot_ui("myModule1"), color_output_box("myModule1")),type = "pills"
          
        )), 
        tabPanel("KEGG Gather", tabsetPanel(tabPanel("Subset Data", subset_data_ui("myModule1")), 
                                            tabPanel("Get Enzyme and Gene Data", KEGG_gather_ui("myModule1")),
                                            tabPanel("Get Pathway Images", pathway_gather_UI("myModule1"))))
      )), 
    mainPanel(tabsetPanel(
      tabPanel("Plot Viewer",tabsetPanel(
        tabPanel("Histogram", verticalLayout(densityOutput("myModule1"),histogramOutput("myModule1"), k_s_table_output("myModule1"))),
        tabPanel("PCA Plot", pcaOutput("myModule1")),
        tabPanel("Random Forest Plots", verticalLayout(varimpOutput("myModule1"), rfPCAOutput("myModule1"))),
        tabPanel("Volcano Plot", verticalLayout(plot_output("myModule1"), table_output("myModule1"))),
        tabPanel("Dotplot", dotplot_output("myModule1")),
        tabPanel("Boxplot", boxplot_output("myModule1")),type = "pills"
        
        
        
      )),
      tabPanel("Data Viewer", tabsetPanel(
        tabPanel("Metabolomics Data",verticalLayout(data_viewer_output("myModule1"), sample_viewer_output("myModule1"))),
        tabPanel("Statistics Data", stats_viewer_output("myModule1")),
        tabPanel("KEGG Gather Data", verticalLayout(KG_viewer_output("myModule1"), genes_viewer_output("myModule1"))))),
      tabPanel("Pathway Viewer", pathway_image_output("myModule1"))
    )
    , width = 6))
)

server <- function(input, output, session) {
  
  
  
  observe({
    # Make sure theme is kept current with desired
    session$setCurrentTheme(
      bslib::bs_theme_update(my_theme, bootswatch = input$current_theme)
    )
  })

  
  omu_list <<- upload_and_normalization_server("myModule1")
  
  PCA_server("myModule1")
  
  random_forest_server("myModule1")
  
  anova_server("myModule1")
  
  volcano_server("myModule1")
  
  dotplot_server("myModule1")
  
  boxplot_server("myModule1")
  
  KEGG_gather_server("myModule1")
  
}

ggplot2::theme_set(ggplot2::theme_bw())
#thematic_shiny()
shinyApp(ui, server)

}

myApp()
