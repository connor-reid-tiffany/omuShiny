#'Function for dynamic color palettes
#' @param id module n
#' @importFrom grDevices hcl

gg_fill_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#' UI for volcano plots
#' @param id module id
#' @importFrom shiny NS actionButton HTML selectizeInput selectInput sliderInput radioButtons downloadButton tagList textInput splitLayout numericInput checkboxInput
#' @importFrom colourpicker colourInput
#' @importFrom shinyFeedback useShinyFeedback
volcano_ui <- function(id){
  
  ns <- NS(id)
  #help button UI
  tagList(
    actionButton(ns("help"), label = "Help", icon = icon("question-circle", lib = "font-awesome",style="color: #fff; background-color: #0694bf; border-color: #013747")),
    HTML("<br><br><h4>Create a Volcano Plot</h4><br><br>"),
    selectInput(ns("stats_data"), "Plot Data", choices = NULL, multiple = FALSE),
    selectInput(ns("fill_variable"), "Metadata Category", choices = c("Class", "Subclass_1", "Subclass_2", "Subclass_3", "Subclass_4")),
    #choices are selected server side based off observed event of fill_variable
    selectizeInput(ns("fill_levels"), "Select Metabolites", choices = NULL,multiple = TRUE),
    #could reduce this to a function, might be useful for other modules
    splitLayout(cellWidths = c("38%", "31%", "31%"),actionButton(ns("exclude_Others"), "Plot Selected",
                                 style="color: #fff; background-color: #0694bf; border-color: #013747"),
                actionButton(ns("include_Others"), "Plot All",
                                 style="color: #fff; background-color: #0694bf; border-color: #013747"),
                actionButton(ns("rev_x"), "Invert FC",style="color: #fff; background-color: #0694bf; border-color: #013747")),
    splitLayout(cellWidths = c("40%", "30%", "30%"),
                numericInput(ns("pval"), "pvalue", value = 0.05,width = "75%", min = 0, max = 1),
                numericInput(ns("l2fc_minus"), "l2fc -", value = -3,width = "70%"),
                numericInput(ns("l2fc_plus"), "l2fc +", value = 3,width = "70%")),
    checkboxInput(ns("label_points"), "label",value = FALSE),
    #could reduce to a function using a dataframe of values for each input and a functional. sliders for plot dimensions etc.
    splitLayout(sliderInput(ns("height"), "Plot Height", min = 100, max = 1500, value = 500, width = "95%"),
                sliderInput(ns("width"), "Plot Width", min = 100, max = 1500, value = 500, width = "95%")),
    splitLayout(numericInput(ns("size"), "Point Size", value = 2, width = "70%"),
                numericInput(ns("font"), "Font Size", value = 9, width = "70%"),
                numericInput(ns("border_size_volcano"), "Border Size", value = 1.5, width = "70%")),
    radioButtons(ns("extension"), "Save As:",
                 choices = c("png", "pdf", "svg", "eps","pptx"), inline = TRUE),
    splitLayout( numericInput(ns("figure_height_volc"), label = "Figure Height(cm)", value = 5, width = "70%"),
                 numericInput(ns("figure_width_volc"), label = "Figure Width(cm)", value = 5, width = "70%")),
    downloadButton(ns("download_plot"), "Save Plot",style="color: #fff; background-color: #0694bf; border-color: #013747")
    
  )
  
}
#'Server for volcano plots
#' @param id module id
#' @importFrom shiny NS moduleServer observeEvent reactiveValues reactiveValuesToList updateSelectizeInput req updateSelectInput reactive renderUI downloadHandler showModal modalDialog renderPlot brushedPoints
#' @importFrom thematic thematic_shiny
#' @importFrom spsComps shinyCatch
#' @importFrom colourpicker colourInput
#' @importFrom ggplot2 ggplot ggsave aes geom_point theme theme_bw element_blank element_text element_rect geom_hline scale_color_manual coord_cartesian 
#' @importFrom ggrepel geom_label_repel
#' @importFrom officer read_pptx ph_with add_slide ph_location_type
volcano_server <- function(id){
  
  moduleServer(id, function(input, output, session){
    thematic::thematic_shiny()
    
    #a reactive to update the data UI based on which data type  is selected as an input
    data_choice <- reactive({
      
      data <- reactiveValuesToList(omu_list)
      stats_data <- data$stats
      stats_data <- stats_data[!names(stats_data)=="Residuals"]
      
      return(stats_data)
      
    })
    
    observeEvent(data_choice(), {
      
      updateSelectInput(inputId = "stats_data", choices = names(data_choice()))
      
    })
    
    data <- reactive({
      req(data_choice())
      
      df <- data_choice()[[input$stats_data]]

      
      return(df)
    })
    
    #a reactive to update the fill_levels UI based on which fill_variable choice is selected as an input
    fill_choice <- reactive({
      fill_data <- data()
      fill_variable <- fill_data[,input$fill_variable]
      return(fill_variable)
      
    })
    
    #an observation which dynamically updates the fill_levels input options when the fill_variable input is changed
    observeEvent(fill_choice(), {
      choices <- unique(fill_choice())
      updateSelectInput(inputId = "fill_levels", choices = c(choices, "Other")) 
    })
    
    #a reactive to make sure the table output from brush is correct
    calc_yaxis <- reactive({
      if (is.null(dat$df)) return()
      calc_data <- dat$df
      calc_data$`-log10(padj)` <- -log10(calc_data$padj)
      calc_data <- calc_data[,!names(calc_data) %in% c("KEGG", "p.value", "null.value","sumsq", 
                                                       "df", "meansq", "estimate", "statistic", 
                                                       "conf.low","conf.high")]
      return(calc_data)
      
    })
    
    #a reactive value meant to change the plot output depending on which create a plot action button is pressed
    #this works by either exluding data that is not in the fill_levels input, or changing exluded data to "Other" to include it
    dat <- reactiveValues(df = NULL)
    
    
    observeEvent(input$include_Others, {
      
      dat$df <- data()
      dat$df[,input$fill_variable] <- sapply(dat$df[, input$fill_variable], function(x) replace(x, !x %in% input$fill_levels, "Other"))
      choices <- unique(fill_choice())
      updateSelectInput(inputId = "fill_levels", choices = c(choices, "Other"), selected = c(input$fill_levels, "Other"))
      data_list <- omu_list$data
      base_L2FC <- calc_base_FCs(base_metabo = data_list$base_metabo, meta = data_list$meta, data_stats = dat$df)
      
      dat$df <- merge(dat$df, base_L2FC, by = "Metabolite")
      
    })
    
    observeEvent(input$exclude_Others, {
      
      dat$df <- data()
      dat$df[,input$fill_variable] <- sapply(dat$df[, input$fill_variable], function(x) replace(x, !x %in% input$fill_levels, NA))
      dat$df <- dat$df[!is.na(dat$df[,input$fill_variable]),]
      data_list <- omu_list$data
      base_L2FC <- calc_base_FCs(base_metabo = data_list$base_metabo, meta = data_list$meta, data_stats = dat$df)
      
      dat$df <- merge(dat$df, base_L2FC, by = "Metabolite")
      
    })
    
    #observation for the invert fold change action button, for if someone wants to swap numerator and denominator contrasts
    observeEvent(input$rev_x, {
      
      if (is.null(dat$df)) return()
      dat$df[,"base_log2FoldChange"] <- dat$df[,"base_log2FoldChange"] * -1
      
      
    })
    
    #reactive value for plot zoom
    ranges <- reactiveValues(x = NULL, y = NULL)
    
    #dynamic UI that only appears once fill_level inputs are selected
    output$myPanel <- renderUI({ 
      lev <- sort(unique(input$fill_levels)) # sorting so that "things" are unambigious
      cols <- gg_fill_hue(length(lev))
      
      # New IDs "colX1" so that it partly coincide with input$select...
      lapply(seq_along(lev), function(i) {
        colourpicker::colourInput(inputId = NS(id,paste0("col", lev[i])),
                                  label = paste0("Choose Color for ", lev[i]), 
                                  value = cbind(cols[i])
        )
        
      })
    })
    
    #reactive which creates the volcano plot. is dependent on data(), and also the reactive value dat(aka plot button being pressed)
    #takes inputs from fill_levels color, font, width, height thus requires these as well
    volcano_output <- reactive(
      {
        if (is.null(dat$df)) return()
        cols <- paste0("c(", paste0("input$col", sort(input$fill_levels), collapse = ", "), ")")
        # print(cols)
        cols <- eval(parse(text = cols))
        # print(cols)
        # To prevent errors
        req(input$fill_levels)
        req(length(cols) == length(input$fill_levels))
        req(data())
        
        #req(input$exclude_Others | input$include_Others
        create_point_labels <- function(stats_data, pval = 1, l2fc_minus= 0, l2fc_plus=0){
          
          labels <- stats_data[stats_data$padj <=pval,]
          labelsminus <- labels[labels$base_log2FoldChange < l2fc_minus,]
          labelsplus <- labels[labels$base_log2FoldChange > l2fc_plus,]
          labels <- rbind(labelsminus, labelsplus)
          labels <- labels[labels$KEGG!="NA",]
          
          return(labels)
          
        }
        
        labels <- create_point_labels(stats_data = dat$df, pval = input$pval, 
                                      l2fc_minus = input$l2fc_minus,l2fc_plus = input$l2fc_plus)
        
        variable <- dat$df[,input$fill_variable]
        if(input$label_points==FALSE){
        ggplot(dat$df, aes(x = base_log2FoldChange, y = -log10(padj), color = unlist(variable), label = Metabolite
                           #,color = unlist(variable), 
                           #shape = unlist(variable)
        )) + 
          geom_point(size = input$size) +
          geom_hline(aes(yintercept = -log10(0.05)),linetype = "dashed", color = "grey", size = 1) +
          #xlim((0-abs(max(dat$df[,"base_log2FoldChange"]))),(0+max(abs(dat$df[,"base_log2FoldChange"]))))+
          scale_color_manual(values = cols) +
          coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = TRUE) +
          #scale_color_manual(values = c(rep("black", length(cols))))+
          #scale_shape_manual(values = c(rep(21, length(cols))))+
          
          theme(panel.grid = element_blank(), panel.border = element_rect(size = input$border_size_volcano),legend.title = element_blank()) +
          theme(axis.text = element_text(size = input$font)) +
          theme(axis.title = element_text(size = input$font)) +
          theme(legend.text = element_text(size = input$font), legend.position = "top")
        }else if (input$label_points==TRUE){
          ggplot(dat$df, aes(x = base_log2FoldChange, y = -log10(padj), colour = unlist(variable), label = Metabolite
                             #,color = unlist(variable), 
                             #shape = unlist(variable)
          )) + 
            geom_point(size = input$size) +
            geom_hline(aes(yintercept = -log10(0.05)),linetype = "dashed", color = "grey", size = 1) +
            #xlim((0-abs(max(dat$df[,"base_log2FoldChange"]))),(0+max(abs(dat$df[,"base_log2FoldChange"]))))+
            geom_label_repel(data = labels,force = 20, colour = "black",
                             fill = "white", show.legend = FALSE,
                             min.segment.length = 0.2, size =4) +
            scale_color_manual(values = cols) +
            coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = TRUE) +
            #scale_color_manual(values = c(rep("black", length(cols))))+
            #scale_shape_manual(values = c(rep(21, length(cols))))+ +
            theme(panel.grid = element_blank(), panel.border = element_rect(size = input$border_size_volcano),legend.title = element_blank()) +
            theme(axis.text = element_text(size = input$font)) +
            theme(axis.title = element_text(size = input$font)) +
            theme(legend.text = element_text(size = input$font), legend.position = "top")
          
        }
      }
    )
    
    #here so that height and width sliders function properly 
    output$volcano <- renderPlot(width = function() input$width, height = function() input$height,res = 96,{
      
      volcano_output()
      
    }, bg = "transparent")
    
    #observer for plot zoom
    observeEvent(input$plot_dblclick, {
      brush <- input$plot_brush
      if (!is.null(brush)) {
        ranges$x <- c(brush$xmin, brush$xmax)
        ranges$y <- c(brush$ymin, brush$ymax)
        
      } else {
        ranges$x <- NULL
        ranges$y <- NULL
      }
    })
    #creates table describing selected datapoints on plot
    output$data <- renderTable({
      req(input$plot_brush)
      brushedPoints(calc_yaxis(), input$plot_brush)
    })
    
    #download handler for plot, the boolean after content is required for pptx to be an option, as it is not covered by ggsave
    output$download_plot <- downloadHandler(
      filename = function() {
        paste("volcano_plot", input$extension, sep = ".")
      },
      
      content = function (file) {
        if (grepl(".pptx", file)==TRUE){
          
          doc <-  read_pptx() 
          doc <- add_slide(doc, layout =  'Title and Content', master = 'Office Theme')
          doc <- ph_with(doc, value = volcano_output(), location = ph_location_type(type = "body")) 
          print(doc, file)
          
        }else {
          volcano_plot <- volcano_output() + theme_bw() + theme(panel.grid = element_blank(), panel.background = element_rect(fill = "white"), 
                               plot.background = element_rect(fill = "white")) + 
            theme(axis.text = element_text(size = input$font),axis.title = element_text(size = input$font),
                  legend.text = element_text(size = input$font),
                  legend.title = element_blank(),
                  panel.grid = element_blank(),
                  panel.border = element_rect(size = input$border_size_volcano), legend.position = "top")
          
          ggsave(file, volcano_plot, device = input$extension,width = input$figure_width_volc, height = input$figure_height_volc)
          
        }
        
      }
    )
    
    #help button pop up text box
    observeEvent(input$help,{
      showModal(modalDialog(
        title = "Help",
        HTML("This is the omu-Shiny module for creating a Volcano plot from data generated in the analysis module.<br>
       <br>
       Perform the following steps to create a plot.<br> 
       <br> 
       1.Import the data.<br> 
       2.Select the metabolite metadata to color by, i.e. Class, Subclass_1, etc. Metadata becomes more descriptive.
       with each successive subclass.<br> 
       3.Select the types of metabolites you want to color.<br> 
       4.Click either 'Create Plot' button to make the plot. 'Create Plot with Selected Metabolites' plots only the points you selected 
       in step 3. 'Create Plot with All Metabolites' includes all datapoints in the plot.<br> 
       <br> 
       Optional Steps:<br> 
       1. To view the inverse foldchange where the numerator and denominator variables are swappped, click the invert foldchange button.<br> 
       2. You can use sliders to adjust plot width and height, point size, and font size.<br> 
       3. You can drag and drop points on the plot to create a table describing what those points are.<br> 
       4. You can zoom in on a drag and dropped section of the plot by double clicking, and zoom out by double clicking again.<br> 
       <br> 
       The plot can be downloaded as either a PNG, PDF, SVG, or PPTX file."
        )))
    })
    
    
    
    
  })
}
#'Function for outputting volcano plot
#' @param id module n
#' @importFrom shiny NS brushOpts plotOutput
plot_output <- function(id){
  ns <- NS(id)
  
  plotOutput(ns("volcano"),dblclick = ns("plot_dblclick"), brush = brushOpts(
    id = ns("plot_brush"),
    resetOnNew = TRUE
  ),width = 1000, height = 500)
  
  
  
}
#'Function for volcano table output
#' @param id module n
#' @importFrom shiny NS tableOutput
table_output <- function(id){
  ns <- NS(id)
  tableOutput(ns("data"))
  
}
#'Function for color output
#' @param id module n
#' @importFrom shiny NS uiOutput
color_output <- function(id){
  
  ns <- NS(id)
  uiOutput(ns("myPanel"))
  
}