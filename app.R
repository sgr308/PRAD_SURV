# ===============================
# Interactive Survival Analysis Shiny App
# ===============================

library(shiny)
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
library(bslib)

# -------------------------------
# UI
# -------------------------------
ui <- fluidPage(
  
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  titlePanel(h2("Gene Expression & Survival Analysis in Prostate Cancer", align = "center")),

  sidebarLayout(
    
    sidebarPanel(
      width = 3,
      
      h4("Data Input"),
      radioButtons(
        "data_mode",
        "Choose data source:",
        choices = c("Use preloaded example data" = "pre",
                    "Upload my own data" = "upload"),
        selected = "pre"
      ),
      
      conditionalPanel(
        condition = "input.data_mode == 'upload'",
        fileInput("topgenes", "Upload TopGenes file", accept = c(".txt", ".csv")),
        fileInput("clinical", "Upload Clinical Survival file", accept = c(".txt", ".csv")),
        fileInput("expr", "Upload Normalized Expression Matrix", accept = c(".txt", ".csv"))
      ),
      
      hr(),
      h4("Gene Selection"),
      selectInput("gene_name", "Select DE Gene:", choices = NULL),
      
      hr(),
      downloadButton("download_pdf", "Download Plots (PDF)")
    ),
    
    mainPanel(
      width = 9,
      
      tabsetPanel(
        
        # -------- HOME TAB --------
        tabPanel("Home",
                 tags$br(),
                 h4("This is an app for exploring bulk RNA-Seq data and evaluating the prognostic impact of gene expression in prostate cancer using TCGA datasets."),
                 tags$br(),
                 h3("Workflow"),
                 tags$img(src = "wrk2.jpg", width = "900"),
                 
                 h3("Benefits"),
                 tags$ul(
                   tags$li("No-code RNA-Seq and survival analysis"),
                   tags$li("Fast biomarker discovery"),
                   tags$li("Interactive and reproducible research")
                 ),
                 h3("App developer"),
                 tags$ul(
                   tags$li(p(tags$b("Sagar Patel, Ph.D."),"|", (tags$b(a(href="https://www.linkedin.com/in/sgr308", target="_blank", "LinkedIn"))))),
                   tags$li(p(tags$b("Email: sgr.bnf(at)gmail.com")))
                 )
        ),
        
        # -------- EXPRESSION TAB --------
        tabPanel("Expression Distribution",
                 br(),
                 plotOutput("boxplot_expr", height = "510px"),
                 br(),
                 uiOutput("expr_text")
        ),
        
        # -------- SURVIVAL TAB --------
        tabPanel("Survival Analysis",
                 br(),
                 plotOutput("km_plot", height = "510px"),
                 br(),
                 uiOutput("surv_text")
        )
      )
    )
  )
)

# -------------------------------
# SERVER
# -------------------------------
server <- function(input, output, session) {
  
  # ---- Load data ----
  load_data <- reactive({
    
    if (input$data_mode == "pre") {
      TopGenes <- read.table("TopGenes_HUGO.txt", header = TRUE, sep = "\t", check.names = FALSE)
      clin_df <- read.table("clinical_data_surv.txt", header = TRUE, sep = "\t", check.names = FALSE)
      d_mat <- read.table("v_normsd.txt", header = TRUE, sep = "\t", check.names = FALSE)
    } else {
      req(input$topgenes, input$clinical, input$expr)
      TopGenes <- read.table(input$topgenes$datapath, header = TRUE, sep = "\t", check.names = FALSE)
      clin_df <- read.table(input$clinical$datapath, header = TRUE, sep = "\t", check.names = FALSE)
      d_mat <- read.table(input$expr$datapath, header = TRUE, sep = "\t", check.names = FALSE)
    }
    
    list(TopGenes = TopGenes, clin_df = clin_df, d_mat = d_mat)
  })
  
  observeEvent(load_data(), {
    updateSelectInput(session, "gene_name",
                      choices = load_data()$TopGenes$Gene_name,
                      selected = load_data()$TopGenes$Gene_name[1])
  })
  
  gene_info <- reactive({
    df <- load_data()$TopGenes
    df[df$Gene_name == input$gene_name, ]
  })
  
  # ---- Expression plot ----
  output$boxplot_expr <- renderPlot({
    gene_id <- gene_info()$Gene_id
    d_mat <- load_data()$d_mat
    clin_df <- load_data()$clin_df
    
    expr_diseased <- d_mat[rownames(clin_df), gene_id]
    expr_healthy <- d_mat[setdiff(rownames(d_mat), rownames(clin_df)), gene_id]
    
    boxplot(expr_diseased, expr_healthy,
            names = c("Diseased", "Healthy"),
            col = c("black", "red"),
            main = paste("Expression of", input$gene_name),
            ylab = "Normalized Expression")
  })
  
  output$expr_text <- renderUI({
    p("Comparison of gene expression between tumor and normal samples.")
  })
  
  # ---- Survival analysis ----
  surv_results <- reactive({
    gene_id <- gene_info()$Gene_id
    clin_df <- load_data()$clin_df
    d_mat <- load_data()$d_mat
    
    clin_df$gene_value <- d_mat[rownames(clin_df), gene_id]
    median_value <- median(clin_df$gene_value, na.rm = TRUE)
    clin_df$gene_group <- ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")
    
    fit <- survfit(Surv(overall_survival, deceased) ~ gene_group, data = clin_df)
    pval <- surv_pvalue(fit, data = clin_df)$pval
    
    list(fit = fit, pval = pval, clin_df = clin_df)
  })
  
  output$km_plot <- renderPlot({
    ggsurvplot(
      surv_results()$fit,
      data = surv_results()$clin_df,
      pval = TRUE,
      risk.table = TRUE,
      title = paste("Survival Analysis:", input$gene_name),
      legend.labs = c("Down-regulated", "Up-regulated")
    )
  })
  
  output$surv_text <- renderUI({
    pval <- surv_results()$pval
    
    if (pval < 0.05) {
      p(strong("Result:"), paste("Significant survival difference (p =", round(pval, 4), ")"))
    } else {
      p(strong("Result:"), paste("No significant survival difference (p =", round(pval, 4), ")"))
    }
  })
  
  # ---- Download PDF ----
  output$download_pdf <- downloadHandler(
    filename = function() paste0(input$gene_name, "_plots.pdf"),
    content = function(file) {
      pdf(file, width = 10, height = 8)
      
      gene_id <- gene_info()$Gene_id
      d_mat <- load_data()$d_mat
      clin_df <- load_data()$clin_df
      
      expr_diseased <- d_mat[rownames(clin_df), gene_id]
      expr_healthy <- d_mat[setdiff(rownames(d_mat), rownames(clin_df)), gene_id]
      
      boxplot(expr_diseased, expr_healthy,
              names = c("Diseased", "Healthy"),
              col = c("black", "red"),
              main = paste("Expression of", input$gene_name))
      
      print(
        ggsurvplot(
          surv_results()$fit,
          data = surv_results()$clin_df,
          pval = TRUE,
          risk.table = TRUE
        )
      )
      
      dev.off()
    }
  )
}

# -------------------------------
# Run app
# -------------------------------
shinyApp(ui = ui, server = server)
