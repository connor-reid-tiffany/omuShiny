gg_fill_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

dotplot_ui <- function(id){
  
  ns <- NS(id)
  #help button UI
  tagList(
    actionButton(ns("help_dotplot"), label = "Help", icon = icon("question-circle", lib = "font-awesome"),style="color: #fff; background-color: #0694bf; border-color: #013747"),
    HTML("<br><br><h4>Create Gene Dotplots</h4><br><br>"),
    HTML("<br><br><h5>Choose Data to Plot</h5><br><br>"),
    selectInput(ns("stats_data_dotplot"), "Select Data to Plot", choices = NULL, multiple = FALSE),
    selectizeInput(ns("Metabolite"), "Select Metabolite to Plot", choices = NULL,multiple = TRUE),
    actionButton(ns("create_plot"), "Create Plot",style="color: #fff; background-color: #0694bf; border-color: #013747"),
    #could reduce to a function using a dataframe of values for each input and a functional. sliders for plot dimensions etc.
    HTML("<br><br><h5>Adjust Plot Parameters</h5><br><br>"),
    splitLayout(sliderInput(ns("height_dotplot"), "Plot Height", min = 100, max = 1500, value = 500),
                sliderInput(ns("width_dotplot"), "Plot Width", min = 100, max = 1500, value = 500)
    ),
    splitLayout(numericInput(ns("size_dotplot"), "Point Size", value = 2, step = 0.25),
                numericInput(ns("font_size_dotplot"), "Font Size", value = 14),
                numericInput(ns("border_size_dotplot"), "Border Size", value = 1.5)
    ),
    colourInput(ns("pvalue_color"), "Choose color for pvalue", value = "black"),
    HTML("<br><br><h5>Download Plot</h5><br><br>"),
    radioButtons(ns("extension_dotplot"), "Save As:",
                 choices = c("pdf", "png","svg", "pptx"), inline = TRUE),
    splitLayout(numericInput(ns("fig_width_dotplot"), "Base Figure Width", value = 5, min = 1, max = 30),
                numericInput(ns("fig_height_dotplot"), "Base Figure Height", value = 5, min = 1, max = 30)
    ),
    downloadButton(ns("download_plot_dotplot"), "Save Plot",style="color: #fff; background-color: #0694bf; border-color: #013747"),
    
    
  )
  
}

dotplot_server <- function(id){
  
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
      
      updateSelectInput(inputId = "stats_data_dotplot", choices = names(data_choice()))
      
    })
    
    data <- reactive({
      req(data_choice())
      
      df <- data_choice()[[input$stats_data_dotplot]]
      df2 <- get_base_data(base_metabo = omu_list$data$base_metabo, meta = omu_list$data$meta, data_stats = df)
      df2$p <- df$padj[match(df2$Metabolite, df$Metabolite)]
      return(df2)
    })
    

    observeEvent(data(), {
      
      
      
      updateSelectInput(inputId = "Metabolite", choices = unique(data()$Metabolite))
      
      
    })
    
    #dynamic UI that only appears once fill_level inputs are selected
    output$myPanel_dotplot <- renderUI({
      lev <- c(unique(data()$term)) # sorting so that "things" are unambigious
      cols <- gg_fill_hue(length(lev))
      
      # New IDs "colX1" so that it partly coincide with input$select...
      lapply(seq_along(lev), function(i) {
        colourpicker::colourInput(inputId = NS(id,paste0("col", lev[i])),
                                  label = paste0("Choose Color for ", lev[i]),
                                  value = cbind(cols[i])
        )
        
      })
    })
    
    dat <- reactiveValues(df = NULL)
    
    
    observeEvent(input$create_plot, {
      
      dat$df <- data()
      dat$df <- dat$df[dat$df$Metabolite %in% input$Metabolite,]
      
      
    })
    
    dot_plot <- reactive({
      
      req(data())
      
      df <- dat$df
      
      # print(cols)
      
      cols <- paste0("c(", paste0("input$col", c(unique(data()$term)), collapse = ", "), ")")
      cols <- eval(parse(text = cols))
      
      pval_df <- data.frame(Metabolite = unique(df$Metabolite), y.position=max(df$Abundance) + 10000000, group1 = levels(as.factor(df$term))[1], 
                            group2 = levels(as.factor(df$term))[2])
      
      pval_df$p <- formatC(df$p[match(pval_df$Metabolite, df$Metabolite)], digits = 2, format = "e")

      plot <- ggplot(data = df, aes(x = term, y = Abundance, fill = term)) + 
        facet_wrap(.~ Metabolite) + 
        geom_jitter(position = position_jitterdodge(), size = input$size_dotplot, shape = 21) + 
        stat_summary(fun = mean, size = 1) + 
        stat_summary(fun.data = mean_se, geom = "errorbar", size = 0.9, width = 0.3)+ 
        coord_trans(y = "log10") + 
        theme(panel.grid = element_blank(), legend.position = "none", panel.border = element_rect(size = input$border_size_dotplot), 
              axis.text = element_text(size = input$font_size_dotplot), strip.text = element_text(size = input$font_size_dotplot), strip.background = element_blank(), 
              axis.title = element_text(size = input$font_size_dotplot), axis.title.x = element_blank()) + scale_fill_manual(values = cols)+
        scale_y_continuous(limits = c(10,1e+08),breaks = c(0.1, 1,10,100,1000,10000,100000,1000000,10000000,100000000), 
                           labels = c(0.1,1, 10,100,1000,10000,100000,1000000,10000000,100000000))
      
       plot <- plot + stat_pvalue_manual(pval_df, label = "p", tip.length = 0.5, size = 4, bracket.size = 0.9,inherit.aes = FALSE, color = input$pvalue_color)
      

        return(plot)

      
    })
    
    dot_plot_save <- reactive({      
      
      req(data())
      
      df <- dat$df
      cols <- paste0("c(", paste0("input$col", sort(unique(data()$term)), collapse = ", "), ")")
      # print(cols)
      cols <- eval(parse(text = cols))
      
      
      pval_df <- data.frame(Metabolite = unique(df$Metabolite), y.position=max(df$Abundance) + 1000000, group1 = levels(as.factor(df$term))[1], 
                            group2 = levels(as.factor(df$term))[2])
      
      pval_df$p <- formatC(df$p[match(pval_df$Metabolite, df$Metabolite)], digits = 2, format = "e")
      
      
      
      plot <- ggplot(data = df, aes(x = term, y = Abundance, fill = term)) + 
        facet_wrap(.~ Metabolite) + 
        geom_jitter(position = position_jitterdodge(), size = input$size_dotplot, shape = 21, color = "black") + 
        stat_summary(fun = mean, size = 1, color = "black") + 
        stat_summary(fun.data = mean_se, geom = "errorbar", size = 0.9, width = 0.3, color = "black")+ 
        coord_trans(y = "log10") +
        theme_bw() +
        theme(panel.grid = element_blank(), panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"), legend.position = "none", panel.border = element_rect(size = input$border_size_dotplot), 
              axis.text = element_text(size = input$font_size_dotplot), strip.text = element_text(size = input$font_size_dotplot), strip.background = element_blank(), 
              axis.title = element_text(size = input$font_size_dotplot), axis.title.x = element_blank()) + scale_fill_manual(values = c(cols, "black"))+
        scale_y_continuous(limits = c(10,1e+08),breaks = c(0.1, 1,10,100,1000,10000,100000,1000000,10000000,100000000), 
                           labels = c(0.1,1, 10,100,1000,10000,100000,1000000,10000000,100000000)) +
        stat_pvalue_manual(pval_df, label = "p", tip.length = 0.5, size = 4, bracket.size = 0.9, inherit.aes = FALSE, color = input$pvalue_color)
      
      
      return(plot)})
    
    
    output$dotplot <- renderPlot(width = function() input$width_dotplot, height = function() input$height_dotplot,res = 96,{
      
      dot_plot()
      
    },bg="transparent")
    
    #creates table describing selected datapoints on plot
    output$dotplot_table <- DT::renderDataTable({
      req(data())
      
      df <- data()
      
      DT::datatable(data = df, options = list(scrollX = TRUE), width = 1800)
    })
    
    output$download_plot_dotplot <- downloadHandler(
      filename = function() {
        paste("dotplot", input$extension_dotplot, sep = ".")
      },
      
      content = function (file) {
        if (grepl(".pptx", file)==TRUE){
          
          doc <-  officer::read_pptx()
          doc <- officer::add_slide(doc, layout =  'Title and Content', master = 'Office Theme')
          doc <- officer::ph_with(doc, value = dot_plot(), location = officer::ph_location_type(type = "body"))
          print(doc, file)
          
        }else {
          
          
          ggsave(file, dot_plot_save(), device = input$extension_dotplot,width = input$fig_width_dotplot, height = input$fig_height_dotplot)
          
        }
        
      }
    )
    
    #help button pop up text box
    observeEvent(input$help,{
      showModal(modalDialog(
        title = "Help",
        HTML("The plot
        on the left shows log2 transformed, median of ratio normalized read counts for each gene by sample, with bars denoting group means and an asterisk
        denoting a signficant difference between groups (if applicable) following FDR correction.
        The plot on the right is a 95% confidence interval of the log2FoldChange between groups, calculated with the following equation:
        log2FoldChange + qnorm(0.025)*lfcSE and log2FoldChange - qnorm(0.025)*lfcSE, for upper and lower bounds respectively. <br>
       <br>
       Perform the following steps to create a plot.<br>
       <br>
        1. Choose criteria to subset genes by. Options are significance, fold change, and KO groups (Optional).<br>
        <br>
        2. Select genes to plot from the drop down menu. Multiple genes can be plotted at once.<br>
        <br>
        The plot can be downloaded as either a PNG, PDF, SVG, or PPTX file."
        ),easyClose = TRUE))
    })
    
    
    
    
    
  })
}

dotplot_output <- function(id){
  ns <- NS(id)
  
  plotOutput(ns("dotplot"))
  
  
  
}

dot_table <- function(id){
  ns <- NS(id)
  DT::dataTableOutput(ns("dotplot_table"))
  
}

color_output_dot <- function(id){
  
  ns <- NS(id)
  uiOutput(ns("myPanel_dotplot"))
  
}