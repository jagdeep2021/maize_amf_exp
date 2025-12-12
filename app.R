# app.R

# ---- Load required packages ----
library(shiny)
library(shinyWidgets)
library(shinyjs)
library(tidyverse)
library(ggplot2)
library(DT)
library(bslib)
library(png)   # for png::readPNG
library(grid)  # for grid::grid.raster
library(zip)   # for zip::zip

# ---- Source the plot functions ----
source("plot_function.R")

# ---- Read the table data once at startup ----
table_data <- read_csv("data/DE_ASE_table.csv") %>%
  select(geneID, Working.Symbol, ASE_type, ASE_sign, Gen_type, Myc_gene)

# ---- UI ----
ui <- fluidPage(
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#2c3e50",
    secondary = "#18bc9c",
    success = "#2ecc71",
    info = "#3498db",
    warning = "#f39c12",
    danger = "#e74c3c"
  ),
  useShinyjs(),
  titlePanel(
    div(
      h1("Maize-AMF Gene Explorer",
         style = "color: #2c3e50; margin-bottom: 20px;"),
      p("Explore gene expression patterns and allele-specific expression in maize with AMF treatment",
        style = "color: #7f8c8d; font-size: 16px;")
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      style = "background-color: #f8f9fa; padding: 20px; border-radius: 5px;",
      div(
        style = "margin-bottom: 20px;",
        h4("Input Parameters", style = "color: #2c3e50; margin-bottom: 15px;"),
        textAreaInput(
          "gene_list",
          "Enter Gene List (one per line):",
          height = "200px",
          placeholder = "Example:\nZm00001eb164860\nbx13"
        ),
        div(
          style = "margin-top: 20px;",
          h5("Filters", style = "color: #2c3e50; margin-bottom: 10px;"),
          checkboxGroupInput(
            "trt_filter", "Select Treatments:",
            choices = c("M", "C"),
            selected = c("M", "C"),
            inline = TRUE
          ),
          checkboxGroupInput(
            "time_filter", "Select Time Points:",
            choices = c("T1", "T2"),
            selected = c("T1", "T2"),
            inline = TRUE
          )
        ),
        div(
          style = "margin-top: 20px;",
          actionButton(
            "plot_button", "Generate Plots",
            class = "btn-primary",
            style = "width: 100%; margin-bottom: 10px;"
          ),
          downloadButton(
            "download_plots", "Download Plots",
            class = "btn-success",
            style = "width: 100%;"
          )
        ),
        progressBar(id = "progress", value = 0, display_pct = TRUE),
        # Simple loading text that shinyjs::show/hide can toggle
        div(
          id = "loading",
          "Generating plots…",
          style = "margin-top: 10px; color: #e74c3c; display: none;"
        )
      )
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel(
          "Results",
          div(
            style = "margin-bottom: 20px;",
            h4("Gene Information", style = "color: #2c3e50;"),
            dataTableOutput("gene_table")
          ),
          div(
            style = "margin-bottom: 20px;",
            h4("Expression Patterns by Genotype", style = "color: #2c3e50;"),
            plotOutput("expression_plot", height = "400px")
          ),
          div(
            style = "margin-bottom: 20px;",
            h4("Allele-Specific Expression Results", style = "color: #2c3e50;"),
            plotOutput("ase_plot", height = "400px")
          )
        )
      )
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {
  # Reactive value to store plot paths
  plots <- reactiveValues(
    de_expr = NULL,
    ase = NULL
  )
  
  # Reactive expression for gene list
  genes <- reactive({
    req(input$gene_list)
    g <- strsplit(input$gene_list, "\n")[[1]]
    g <- g[g != ""]
    g
  })
  
  # Generate plots when button is clicked
  observeEvent(input$plot_button, {
    show("loading")
    updateProgressBar(session, "progress", value = 0,
                      title = "Starting plot generation...")
    
    genes_list <- genes()
    
    if (length(genes_list) == 0) {
      showNotification("Please enter at least one gene", type = "error")
      hide("loading")
      return()
    }
    
    tryCatch({
      # ---- Expression Patterns by Genotype ----
      updateProgressBar(session, "progress", value = 20,
                        title = "Generating expression plot...")
      
      expr_ok <- tryCatch({
        res <- plot_expression_and_amreads(
          query_genes = genes_list,
          title_prefix = "Expression Patterns by Genotype",
          trt_filter = input$trt_filter,
          time_filter = input$time_filter
        )
        # res is a list with expr_plot + amreads_plot
        plots$de_expr <- res$expr_plot
        TRUE
      }, error = function(e) {
        showNotification(
          paste("Error generating expression plot:", e$message),
          type = "error"
        )
        FALSE
      })
      
      if (!expr_ok) {
        updateProgressBar(session, "progress", value = 0,
                          title = "Error in expression plot")
        hide("loading")
        return()
      }
      
      # ---- ASE plot ----
      updateProgressBar(session, "progress", value = 60,
                        title = "Generating ASE plot...")
      
      ase_ok <- tryCatch({
        res2 <- plot_ase_expression(
          query_genes = genes_list,
          title_prefix = "Allele-Specific Expression Results",
          trt_filter = input$trt_filter,
          time_filter = input$time_filter
        )
        # res2 is TRUE in your current function; ASE plot path is deterministic
        plots$ase <- file.path("plots",
                               "Allele-Specific Expression Results_ASE_expression_plot.png")
        TRUE
      }, error = function(e) {
        showNotification(
          paste("Error generating ASE plot:", e$message),
          type = "error"
        )
        FALSE
      })
      
      if (!ase_ok) {
        updateProgressBar(session, "progress", value = 0,
                          title = "Error in ASE plot")
        hide("loading")
        return()
      }
      
      updateProgressBar(session, "progress", value = 100,
                        title = "Done!")
      hide("loading")
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
      updateProgressBar(session, "progress", value = 0,
                        title = "Error occurred")
      hide("loading")
    })
  })
  
  # Render expression plot
  output$expression_plot <- renderPlot({
    req(plots$de_expr)
    tryCatch({
      if (file.exists(plots$de_expr)) {
        img <- png::readPNG(plots$de_expr)
        grid::grid.raster(img)
      } else {
        plot.new()
        text(0.5, 0.5,
             paste("Expression plot file not found:\n", plots$de_expr))
      }
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5,
           paste("Error loading expression plot:\n", e$message))
    })
  })
  
  # Render ASE plot
  output$ase_plot <- renderPlot({
    req(plots$ase)
    tryCatch({
      if (file.exists(plots$ase)) {
        img <- png::readPNG(plots$ase)
        grid::grid.raster(img)
      } else {
        plot.new()
        text(0.5, 0.5,
             paste("ASE plot file not found:\n", plots$ase))
      }
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5,
           paste("Error loading ASE plot:\n", e$message))
    })
  })
  
  # Gene info table
  output$gene_table <- renderDataTable({
    selected_genes <- unlist(strsplit(input$gene_list, "\n"))
    selected_genes <- selected_genes[selected_genes != ""]
    
    if (length(selected_genes) == 0) {
      return(NULL)
    }
    
    table_data %>%
      filter(geneID %in% selected_genes | Working.Symbol %in% selected_genes)
  }, options = list(
    pageLength = 10,
    scrollX = TRUE,
    autoWidth = TRUE,
    searching = TRUE
  ))
  
  # Download handler for plots
  output$download_plots <- downloadHandler(
    filename = function() {
      "plots.zip"
    },
    content = function(file) {
      tryCatch({
        plot_files <- c(plots$de_expr, plots$ase)
        plot_files <- plot_files[!is.null(plot_files)]
        
        if (length(plot_files) == 0 || !all(file.exists(plot_files))) {
          stop("Some plot files are missing. Please generate the plots first.")
        }
        
        temp_dir <- tempdir()
        file.copy(plot_files, temp_dir, overwrite = TRUE)
        
        zip_file <- file.path(temp_dir, "plots.zip")
        zip::zip(zip_file,
                 files = basename(plot_files),
                 root = temp_dir)
        
        file.copy(zip_file, file)
      }, error = function(e) {
        showNotification(
          paste("Error creating zip file:", e$message),
          type = "error"
        )
      })
    }
  )
  
  # Clear plot paths when inputs change
  observe({
    input$gene_list
    input$trt_filter
    input$time_filter
    plots$de_expr <- NULL
    plots$ase <- NULL
  })
}

# ---- Run app ----
shinyApp(ui = ui, server = server)