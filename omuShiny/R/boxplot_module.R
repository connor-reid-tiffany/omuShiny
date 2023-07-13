gg_fill_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


boxplot_ui <- function(id){
  
  ns <- NS(id)
  #help button UI
  tagList(
    actionButton(ns("help_boxplot"), label = "Help", icon = icon("question-circle", 
    lib = "font-awesome"),style="color: #fff; background-color: #0694bf; border-color: #013747"),
    HTML("<br><br><h4>Create Metabolite Boxplots</h4><br><br>"),
    HTML("<br><br><h5>Choose Data to Plot</h5><br><br>"),
    selectInput(ns("stats_data_boxplot"), "Select Data to Plot", choices = NULL, multiple = FALSE),
    selectizeInput(ns("Metabolite_boxplot"), "Select Metabolite to Plot", choices = NULL,multiple = TRUE),
    actionButton(ns("create_plot_boxplot"), "Create Plot",
    style="color: #fff; background-color: #0694bf; border-color: #013747"),
    #could reduce to a function using a dataframe of values for each input and a functional. sliders for plot dimensions etc.
    HTML("<br><br><h5>Adjust Plot Parameters</h5><br><br>"),
    splitLayout(sliderInput(ns("height_boxplot"), "Plot Height", min = 100, max = 1500, value = 500),
                sliderInput(ns("width_boxplot"), "Plot Width", min = 100, max = 1500, value = 500)
    ),
    splitLayout(numericInput(ns("size_boxplot"), "Point Size", value = 2, step = 0.25),
                numericInput(ns("font_size_boxplot"), "Font Size", value = 14),
                numericInput(ns("border_size_boxplot"), "Border Size", value = 1.5)
    ),
    colourInput(ns("pvalue_color_box"), "Choose color for pvalue", value = "black"),
    HTML("<br><br><h5>Download Plot</h5><br><br>"),
    radioButtons(ns("extension_boxplot"), "Save As:",
                 choices = c("pdf", "png","svg","eps", "pptx"), inline = TRUE),
    splitLayout(numericInput(ns("fig_width_boxplot"), "Base Figure Width", value = 5, min = 1, max = 30),
                numericInput(ns("fig_height_boxplot"), "Base Figure Height", value = 5, min = 1, max = 30)
    ),
    downloadButton(ns("download_plot_boxplot"), "Save Plot",
                   style="color: #fff; background-color: #0694bf; border-color: #013747"),
    
    
  )
  
}

boxplot_server <- function(id){
  
  moduleServer(id, function(input, output, session){
    thematic::thematic_shiny()
    
    
    #a reactive to update the data UI based on which data type  is selected as an input
    data_choice_box <- reactive({
      
      data <- reactiveValuesToList(omu_list)
      stats_data <- data$stats
      stats_data <- stats_data[!names(stats_data)=="Residuals"]
      
      return(stats_data)
      
    })
    
    observeEvent(data_choice_box(), {
      
      updateSelectInput(inputId = "stats_data_boxplot", choices = names(data_choice_box()))
      
    })
    
    data_box <- reactive({
      req(data_choice_box())
      
      df <- data_choice_box()[[input$stats_data_boxplot]]
      df2 <- get_base_data(base_metabo = omu_list$data$base_metabo, meta = omu_list$data$meta, data_stats = df)
      df2$p <- df$padj[match(df2$Metabolite, df$Metabolite)]
      return(df2)
    })
    
    
    observeEvent(data_box(), {
      
      
      
      updateSelectInput(inputId = "Metabolite_boxplot", choices = unique(data_box()$Metabolite))
      
      
    })
    
    #dynamic UI that only appears once fill_level inputs are selected
    output$myPanel_boxplot <- renderUI({
      lev <- c(unique(data_box()$term)) # sorting so that "things" are unambigious
      cols <- gg_fill_hue(length(lev))
      
      # New IDs "colX1" so that it partly coincide with input$select...
      lapply(seq_along(lev), function(i) {
        colourpicker::colourInput(inputId = NS(id,paste0("col_box", lev[i])),
                                  label = paste0("Choose Color for ", lev[i]),
                                  value = cbind(cols[i])
        )
        
      })
    })
    
    dat_box <- reactiveValues(df = NULL)
    
    
    observeEvent(input$create_plot_boxplot, {
      
      dat_box$df <- data_box()
      dat_box$df <- dat_box$df[dat_box$df$Metabolite %in% input$Metabolite_boxplot,]
      
      
    })
    
    #function to generate breaks
    every_nth = function(n) {
      return(function(x) {x[c(TRUE, rep(FALSE, n - 1))]})
    }
    
    box_plot <- reactive({
      
      req(data_box())
      
      df <- dat_box$df
      
      # print(cols)
      
      cols <- paste0("c(", paste0("input$col_box", c(unique(data_box()$term)), collapse = ", "), ")")
      cols <- eval(parse(text = cols))
      
      pval_df <- data.frame(Metabolite = unique(df$Metabolite), y.position=max(df$Abundance) - 50000, group1 = levels(as.factor(df$term))[1], 
                            group2 = levels(as.factor(df$term))[2])
      df_list <- split(df, f = as.factor(df$Metabolite))
      
      max_abundance_list <- lapply(df_list, function(x) max(x$Abundance))
      
      max_abundance_df <- as.data.frame(do.call("rbind", max_abundance_list))
      
      colnames(max_abundance_df)[1] <- "Abundance" 
      
      max_abundance_df$Metabolite <- rownames(max_abundance_df)
      
      pval_df$y.position <- max_abundance_df$Abundance[match(pval_df$Metabolite, max_abundance_df$Metabolite)] 
      
      pval_df$y.position <- pval_df$y.position - (0.1 * pval_df$y.position)
      
      pval_df$p <- formatC(df$p[match(pval_df$Metabolite, df$Metabolite)], digits = 2, format = "e")
      
      plot <- ggplot(data = df, aes(x = term, y = Abundance, fill = term)) + 
        facet_wrap(.~ Metabolite, scales = "free") + 
        geom_boxplot(alpha = 1/3, size = input$size_boxplot) +
        scale_fill_manual(values = cols) +
        theme(panel.grid = element_blank(), legend.position = "none", 
              panel.border = element_rect(size = input$border_size_boxplot), 
              axis.text = element_text(size = input$font_size_boxplot), 
              strip.text = element_text(size = input$font_size_boxplot), 
              strip.background = element_blank(), 
              axis.title = element_text(size = input$font_size_boxplot), axis.title.x = element_blank())
      
      plot <- plot + stat_pvalue_manual(pval_df, label = "p", tip.length = 0.02, size = 3.5, 
                    bracket.size = 0.9,inherit.aes = FALSE, color = input$pvalue_color_box, 
                    bracket.shorten = 0.5)
      
      
      return(plot)
      
      
    })
    
    box_plot_save <- reactive({      
      
      req(data_box())
      
      df <- dat_box$df
      cols <- paste0("c(", paste0("input$col_dot", sort(unique(data_box()$term)), collapse = ", "), ")")
      # print(cols)
      cols <- eval(parse(text = cols))
      
      
      pval_df <- data.frame(Metabolite = unique(df$Metabolite), y.position=max(df$Abundance) - 50000, group1 = levels(as.factor(df$term))[1], 
                            group2 = levels(as.factor(df$term))[2])
      
      df_list <- split(df, f = as.factor(df$Metabolite))
      
      max_abundance_list <- lapply(df_list, function(x) max(x$Abundance))
      
      max_abundance_df <- as.data.frame(do.call("rbind", max_abundance_list))
      
      colnames(max_abundance_df)[1] <- "Abundance" 
      
      max_abundance_df$Metabolite <- rownames(max_abundance_df)
      
      pval_df$y.position <- max_abundance_df$Abundance[match(pval_df$Metabolite, max_abundance_df$Metabolite)] 
      
      pval_df$y.position <- pval_df$y.position - (0.1 * pval_df$y.position)
      
      pval_df$p <- formatC(df$p[match(pval_df$Metabolite, df$Metabolite)], digits = 2, format = "e")
      
      
      
      plot <- ggplot(data = df, aes(x = term, y = Abundance, fill = term)) + 
        facet_wrap(.~ Metabolite, scales = "free") + 
        geom_boxplot(size = input$size_boxplot, alpha = 1/3) + 
        scale_fill_manual(values = cols) +
        theme_bw() +
        theme(panel.grid = element_blank(), panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"), legend.position = "none", panel.border = element_rect(size = input$border_size_dotplot), 
              axis.text = element_text(size = input$font_size_boxplot), strip.text = element_text(size = input$font_size_boxplot), strip.background = element_blank(), 
              axis.title = element_text(size = input$font_size_boxplot), axis.title.x = element_blank()) +
        stat_pvalue_manual(pval_df, label = "p", tip.length = 0.02, size = 3.5, bracket.size = 0.9, 
                           inherit.aes = FALSE, color = input$pvalue_color_box, 
                           bracket.shorten = 0.5)
      
      
      return(plot)})
    
    
    output$boxplot <- renderPlot(width = function() input$width_boxplot, height = function() input$height_boxplot,
                                 res = 96,{
      
      box_plot()
      
    },bg="transparent")
    
    #creates table describing selected datapoints on plot
    output$boxplot_table <- DT::renderDataTable({
      req(data_box())
      
      df <- data_box()
      
      DT::datatable(data = df, options = list(scrollX = TRUE), width = 1800)
    })
    
    output$download_plot_boxplot <- downloadHandler(
      filename = function() {
        paste("dotplot", input$extension_boxplot, sep = ".")
      },
      
      content = function (file) {
        if (grepl(".pptx", file)==TRUE){
          
          doc <-  officer::read_pptx()
          doc <- officer::add_slide(doc, layout =  'Title and Content', master = 'Office Theme')
          doc <- officer::ph_with(doc, value = box_plot(), location = officer::ph_location_type(type = "body"))
          print(doc, file)
          
        }else {
          
          
          ggsave(file, box_plot_save(), device = input$extension_boxplot,width = input$fig_width_boxplot, height = input$fig_height_boxplot)
          
        }
        
      }
    )
    
    #help button pop up text box
    observeEvent(input$help,{
      showModal(modalDialog(
        title = "Help",
        HTML("The plot shows boxplots of metabolites with pvalues. 
        Requires data from the univariate statistics module. If you used a parametric test, use the dotplot module instead.<br>
       <br>
       Perform the following steps to create a plot.<br>
       <br>
        1. Choose criteria to subset metabolites by. Options are significance, fold change, and KO groups (Optional).<br>
        <br>
        2. Select metabolites to plot from the drop down menu. Multiple metabolites can be plotted at once.<br>
        <br>
        The plot can be downloaded as either a PNG, PDF, SVG, or PPTX file."
        ),easyClose = TRUE))
    })
    
    
    
    
    
  })
}

boxplot_output <- function(id){
  ns <- NS(id)
  
  plotOutput(ns("boxplot"))
  
  
  
}

box_table <- function(id){
  ns <- NS(id)
  DT::dataTableOutput(ns("boxplot_table"))
  
}

color_output_box <- function(id){
  
  ns <- NS(id)
  uiOutput(ns("myPanel_boxplot"))
  
}