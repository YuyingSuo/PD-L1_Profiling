
library(shiny)
library(shinyWidgets)

library(shinycssloaders)
library(shinyjs)
library(markdown)

app_source_dir <- file.path("R")
source(file.path(app_source_dir, "mod_plot.R"))
source(file.path(app_source_dir, "helper.R"))


err_txt2='Error: The specified Gene does not exist in this data matrix.'
err_txt3='Error: The specified PhosphoSite does not exist in this data matrix.'

data_dir <- file.path("Data")
load(file.path(data_dir, 'KEGG_1.42.0_genes.RData'))
gmt=dt

########################################################################################################
# User Interface     User Interface     User Interface     User Interface     User Interface     
########################################################################################################

ui <- fluidPage(
  
  #@@@@@@@@@@@@@@@@@@@@@
  ### Prepare the Header Panel
  
  tags$div(
    style = "padding: 10px;",
    tags$img(src = "pdl1-2.png", style = "width: 100%; max-width: 400px; height: auto;")
  ),
  hr(),
  
  
  #@@@@@@@@@@@@@@@@@@@@@
  ### Prepare the side Panel
  
  # adjust the width 
  tags$head(
    tags$style(HTML("
    .responsive-sidebar {
        width: 100%;max-width: 300px;min-width: 200px;
      }
    "))
  ),
  
  
  
  sidebarLayout(
    sidebarPanel(
      class = "responsive-sidebar", 
      
      p("Documents for Getting Started:"),
      downloadLink("Manual", label= "Download the User Manual PDF"),
      hr(),
      
      radioButtons("cancer", label = "Cancer type:",
                   choices =list(               # choices : cancertypes
                     "AML","ccRCC","COAD","GBM","HCC","HNSCC","LSCC","LUAD","OV","PDAC","SCLC","UCEC"),
                   selected = 1),
      radioButtons("datatype", label = "Choose the dataset type:", 
                   choices = list('Proteomic','Transcriptomic'),
                   selected = 1),
    ),
    
    
    #@@@@@@@@@@@@@@@@@@@@@
    ### Prepare the Main Data Panel
    
    mainPanel(
      tabsetPanel(type = "tabs",
                  id = "main_tab", 
                  tabPanel("Welcome", includeMarkdown("Welcome.Rmd")),
                  
                  tabPanel("Expression", 
                           br(),
                           uiOutput("exp_error"),
                           withSpinner(plotOutput(outputId = "expr",width = "600px", height = "500px")),
                           br(),
                           downloadButton("expr_plot", "Download_Plot"),
                           downloadButton("expr_data", "Download_Dataset")), 
                  
                  tabPanel("DEGs", 
                           br(),
                           withSpinner(uiOutput("deg_error")),
                           downloadButton("deg_data_all", "Download_All_DEGs"),
                           br(),
                           textInput("gene", "Enter Gene Symbol"),
                           actionButton("deg_go", "Run Test"),
                           uiOutput("deg_error2"),
                           plotOutput(outputId = "deg",width = "600px", height = "500px"),
                           br(),
                           downloadButton("deg_plot", "Download_Plot"),
                           downloadButton("deg_data", "Download_Dataset")),
                  
                  tabPanel("Correlations", 
                           br(),
                           withSpinner(uiOutput("cor_error")),
                           downloadButton("cor_all", "Download_ALL_correlation data_Sig"),
                           
                           textInput("gene2", "Enter Gene Symbol"),
                           actionButton("cor_go", "Run Test"),
                           uiOutput("cor_error2"),
                           plotOutput(outputId = "cor",width = "600px", height = "500px"),
                           br(),
                           downloadButton("cor_plot", "Download_Plot"),
                           downloadButton("cor_data", "Download_Dataset")), 
                  
                  tabPanel("Pathway Activity", 
                           br(),
                           uiOutput("path_error"),
                           withSpinner(plotOutput(outputId = "ssgsea",width = "600px", height = "500px")),
                           br(),
                           downloadButton("ssgsea_plot", "Download_Plot"),
                           downloadButton("ssgsea_data", "Download_Dataset")),
                  
                  tabPanel("Pathway crosstalk",  
                           br(),
                           uiOutput("cs_notice"),
                           withSpinner(downloadButton("path_cs_all", "Download_ALL_Crosstalk_Sig")),
                           br(),
                           selectizeInput("path_cs", "Choose the pathway:", choices =  c("Please select a pathway" = "", names(gmt)),
                                          selected = "", multiple = FALSE),
                           plotOutput(outputId = "cs",width = "600px", height = "500px"),
                           br(),
                           downloadButton("cs_plot", "Download_Plot"),
                           downloadButton("cs_data", "Download_Dataset")), 
                  
                  tabPanel("Phosphorylation", 
                           br(),
                           uiOutput("pho_notice"),
                           downloadButton("pho_sig_all", "Download_ALL_Phosphosites_Sig"),
                           br(),
                           textInput("phosite", "Enter PhosphoSites"),
                           uiOutput("pho_inputnotice"),  
                           actionButton("pho_go", "Run Test"),
                           uiOutput("pho_error"),
                           br(),
                           withSpinner(plotOutput(outputId = "pho",width = "600px", height = "500px")),
                           downloadButton("pho_plot", "Download_Plot"),
                           downloadButton("pho_data", "Download_Dataset")) 
      ),
      tags$div(
        "PD-L1 Multi-Omics Profiling, version 1.0 (2025)",
        style = "text-align: center; margin-top: 500px; color: #666;")
    )
  )
)






########################################################################################################
# Server     Server     Server     Server     Server     Server     Server     Server     Server     
########################################################################################################

server <- function(input, output, session) {
  # === Global heartbeat: Refresh connection every 20 seconds. ===
  observe({
    invalidateLater(30000, session)  
    session$sendCustomMessage("keepalive", list(t = Sys.time()))
  })
  
  #========================
  #  Lightweight cache and cleaning tool
  #========================
  
  # "lightweight result" cache (PNG + small table)
  .expr_result_png  <- reactiveVal(NULL)
  .deg_result_png   <- reactiveVal(NULL)
  .cor_result_png   <- reactiveVal(NULL)
  .path_result_png  <- reactiveVal(NULL)
  .cs_result_png    <- reactiveVal(NULL)
  .pho_result_png   <- reactiveVal(NULL)
  
  .expr_export_df <- reactiveVal(NULL)
  .deg_export_df  <- reactiveVal(NULL)
  .cor_export_df  <- reactiveVal(NULL)
  .path_export_df <- reactiveVal(NULL)
  .cs_export_df   <- reactiveVal(NULL)
  .pho_export_df  <- reactiveVal(NULL)
  
  .path_cs_all_df <- reactiveVal(NULL)
  .pho_sig_all_df <- reactiveVal(NULL)
  .pdl1_vec <- reactiveVal(NULL)
  .ssgsea_t_cache <- reactiveVal(NULL)
  .pho_site_last <- reactiveVal(NULL)
  .pho_full_df    <- reactiveVal(NULL)   
  .pho_sig_all_df <- reactiveVal(NULL)
  
  # Each module cache / flag
  deg_cache <- reactiveVal(list())
  previous_input_deg <- reactiveVal(NULL)
  
  cor_cache <- reactiveVal(list())
  previous_input_cor <- reactiveVal(NULL)
  
  ssgsea_cache <- reactiveVal(NULL)
  previous_input <- reactiveVal(NULL)
  
  # Large object cleanup + Forced GC
  .clear_big <- function(...) {
    for (x in list(...)) if (exists(x, inherits = TRUE)) rm(list = x, inherits = TRUE)
    invisible(gc())
  }
  
  # TTL: If necessary, the lightweight cache can be cleared at regular intervals to prevent idle usage from occupying RAM.
  .schedule_ttl <- function(val_fun, seconds = 300) {
    later::later(function() { val_fun(NULL); invisible(gc()) }, delay = seconds)
  }
  
  # Secure drawing: Write to temporary file → Read in → Delete
  .plot_to_raw <- function(plot_expr) {
    tf <- tempfile(fileext = ".png")
    on.exit({ if (file.exists(tf)) unlink(tf) }, add = TRUE)
    grDevices::png(filename = tf)
    try({
      plot_expr()
    }, silent = TRUE)
    grDevices::dev.off()
    readBin(tf, what = "raw", n = file.size(tf))
  }
  
  # write a data.frame to CSV in row-chunks to avoid huge string buffers
  .write_csv_chunked <- function(df, file, chunk_rows = 20000, row_names = TRUE) {
    n <- nrow(df)
    if (is.null(n) || n == 0) {
      write.csv(data.frame(Message = "No data to export"), file, row.names = FALSE)
      return(invisible())
    }
    # first chunk (write header)
    i <- 1L
    j <- min(n, chunk_rows)
    write.table(df[i:j, , drop = FALSE], file = file, sep = ",",
                row.names = row_names, col.names = NA, qmethod = "double")
    i <- j + 1L
    # append remaining chunks without header
    while (i <= n) {
      j <- min(n, i + chunk_rows - 1L)
      write.table(df[i:j, , drop = FALSE], file = file, sep = ",",
                  row.names = row_names, col.names = FALSE, qmethod = "double",
                  append = TRUE)
      i <- j + 1L
    }
  }
  
  # safer getter for reactiveVal that won't bind reactive deps
  .get_isolated <- function(rv) {
    isolate(rv())
  }
  
  # When switching to the main tab, clear the light cache uniformly.
  observeEvent(input$main_tab, {
    .expr_result_png(NULL); .deg_result_png(NULL); .cor_result_png(NULL)
    .path_result_png(NULL); .cs_result_png(NULL);  .pho_result_png(NULL)
    .expr_export_df(NULL);  .deg_export_df(NULL);  .cor_export_df(NULL)
    .path_export_df(NULL);  .cs_export_df(NULL);   .pho_export_df(NULL)
    .path_cs_all_df(NULL);  .pho_sig_all_df(NULL);  .pdl1_vec <- reactiveVal(NULL)
    .pho_site_last(NULL);  .pho_full_df(NULL);  .pho_sig_all_df(NULL)
    invisible(gc())
  }, ignoreInit = TRUE)
  
  # Input change: Clear all results and memory status, ensuring a fresh analysis.
  observeEvent(list(input$cancer, input$datatype), {
    .expr_result_png(NULL); .deg_result_png(NULL); .cor_result_png(NULL)
    .path_result_png(NULL); .cs_result_png(NULL);  .pho_result_png(NULL)
    .expr_export_df(NULL);  .deg_export_df(NULL);  .cor_export_df(NULL)
    .path_export_df(NULL);  .cs_export_df(NULL);   .pho_export_df(NULL)
    .path_cs_all_df(NULL);  .pho_sig_all_df(NULL);  .pdl1_vec <- reactiveVal(NULL)
    .ssgsea_t_cache <- reactiveVal(NULL); .pho_site_last(NULL)
    .pho_full_df(NULL);  .pho_sig_all_df(NULL)
    deg_cache(list())
    cor_cache(list())
    ssgsea_cache(NULL)
    
    previous_input_deg(NULL)
    previous_input_cor(NULL)
    previous_input(NULL)
    
    invisible(gc())
  }, ignoreInit = TRUE)
  
  # Session cleanup
  session$onSessionEnded(function() {
    .expr_result_png(NULL); .deg_result_png(NULL); .cor_result_png(NULL)
    .path_result_png(NULL); .cs_result_png(NULL);  .pho_result_png(NULL)
    .expr_export_df(NULL);  .deg_export_df(NULL);  .cor_export_df(NULL)
    .path_export_df(NULL);  .cs_export_df(NULL);   .pho_export_df(NULL)
    .path_cs_all_df(NULL);  .pho_sig_all_df(NULL);  .pdl1_vec <- reactiveVal(NULL)
    .ssgsea_t_cache <- reactiveVal(NULL); .pho_site_last(NULL)
    .pho_full_df(NULL);  .pho_sig_all_df(NULL)
    invisible(gc())
  })
  
  # Lightweight periodic GC
  manage_memory <- local({
    rt <- shiny::reactiveTimer(30000)
    force(function() {
      observe({
        rt()
        invisible(gc())
      })
    })
  })
  manage_memory()
  
  
  #========================
  #  user mannual
  #========================
  output$Manual <- downloadHandler(
    filename = "User Manual.pdf",
    content = function(file){
      file.copy("User Manual.pdf", file)
    }
  )
  
  
  #===========================================================
  #  Expression
  #===========================================================
  
  exp_test <- reactive({
    req(input$main_tab == "Expression")
    req(input$cancer, input$datatype)
    
    cn   <- input$cancer
    type <- input$datatype
    
    dt1 <- dataprocess1(cn, type) 
    expr_data <- dt1$exp
    err       <- dt1$err
    
    if (is.null(expr_data)) {
      expr_row <- NA
    } else {
      idx <- which(rownames(expr_data) == "CD274")
      expr_row <- if (length(idx) == 0) NA else expr_data[idx, , drop = FALSE]
    }
    
    .expr_export_df(if (all(is.na(expr_row))) NULL else t(expr_row))
    
    # Pre-rendered screen display/export in PNG bytes (to temporary file)
    if (!all(is.na(expr_row))) {
      bytes <- .plot_to_raw(function() {
        if (type == "Transcriptomic" && (cn %in% rna_noNAT)) {
          print(draw_expr_plot_onlyT(expr_row))
        } else {
          print(draw_expr_plot(expr_row))
        }
      })
      .expr_result_png(bytes)
      .schedule_ttl(.expr_result_png)
    } else {
      .expr_result_png(NULL)
    }
    
    .clear_big("dt1", "expr_data")
    
    return(list(exp = expr_row, err = err))
  })
  
  output$exp_error <- renderUI({
    req(input$main_tab == "Expression")
    res <- exp_test()
    if (!is.null(res$err)) {
      tags$div(style = "color: red; font-weight: bold;", res$err)
    }
  })
  
  output$expr_data <- downloadHandler(
    filename = function() {
      paste0(input$cancer,"_",input$datatype,"_PD-L1_expression_", Sys.Date(), ".csv")
    },
    content = function(file) {
      df2 <- .expr_export_df()
      if (is.null(df2)) {
        write.csv(data.frame(Message = "No data to export"), file, row.names = FALSE)
      } else {
        write.csv(df2, file, row.names = TRUE)
      }
      .expr_export_df(NULL); invisible(gc())
    }
  )
  
  output$expr <- renderPlot({
    req(input$main_tab == "Expression")
    res <- exp_test()
    if (all(is.na(res$exp))) return(NULL)
    cn <- input$cancer; type <- input$datatype
    if (type == "Transcriptomic" && (cn %in% rna_noNAT)) {
      draw_expr_plot_onlyT(res$exp)
    } else {
      draw_expr_plot(res$exp)
    }
  })
  
  output$expr_plot <- downloadHandler(
    filename = function() {
      paste0(input$cancer,"_",input$datatype,"_PD-L1_expression_", Sys.Date(),".png")
    },
    content = function(file) {
      bytes <- .expr_result_png()
      if (is.null(bytes)) {
        res <- exp_test()
        if (all(is.na(res$exp))) return(NULL)
        tf <- tempfile(fileext = ".png")
        on.exit({ if (file.exists(tf)) unlink(tf) }, add = TRUE)
        grDevices::png(tf)
        cn <- input$cancer; type <- input$datatype
        if (type == "Transcriptomic" && (cn %in% rna_noNAT)) {
          print(draw_expr_plot_onlyT(res$exp))
        } else {
          print(draw_expr_plot(res$exp))
        }
        grDevices::dev.off()
        file.copy(tf, file)
      } else {
        writeBin(bytes, file)
      }
      .expr_result_png(NULL); invisible(gc())
    }
  )
  
  
  
  #===========================================================
  #  DEGs
  #===========================================================
  
  observeEvent(list(input$main_tab,input$cancer, input$datatype), {
    req(input$main_tab=="DEGs")
    req(input$cancer,input$datatype)
    
    current_input_deg <- list(a = input$cancer, b = input$datatype)
    test <- is.null(previous_input_deg()) || !identical(current_input_deg, previous_input_deg())
    if (!test) return()
    
    previous_input_deg(current_input_deg)
    cn <- input$cancer; type <- input$datatype
    
    dt2 <- dataprocess1(cn, type)
    expr_data <- dt2$exp
    err <- dt2$err
    
    if (is.null(expr_data)) {
      deg_sig <- NA
      cd274_2group <- NA
      deg_cache(list(cd274_2group=cd274_2group, deg_sig = deg_sig, err = err))
      return(invisible(NULL))
    }
    
    if (type=='Transcriptomic' && (cn %in% rna_noNAT)) {
      expr_t <- expr_data
    } else {
      expr_t <- expr_data[, -grep("_N", colnames(expr_data)), drop = FALSE]
    }
    
    cd274_2group <- deg_dataprocess(expr_t)    
    deg_sig <- cd274_2group[which(cd274_2group$change_padj != 'NoSignifi'), , drop = FALSE]
    deg_cache(list(cd274_2group=cd274_2group, deg_sig = deg_sig, err = NULL))
    
    .clear_big("dt2", "expr_data", "expr_t")
  }, ignoreInit = TRUE)
  
  output$deg_error <- renderUI({
    req(input$main_tab == "DEGs")
    res2 <- deg_cache()
    if (!is.null(res2$err)) {
      tags$div(style = "color: red; font-weight: bold;", res2$err)
    }
  })
  
  
  select_deg <- eventReactive(input$deg_go, {
    gene <- toupper(input$gene)
    
    # isolate cache to avoid recompute/extra copies
    ct <- .get_isolated(deg_cache)
    sigdata <- ct$cd274_2group
    
    if (is.null(sigdata) || is.null(rownames(sigdata))) {
      .deg_result_png(NULL); .deg_export_df(NULL)
      return(list(deg_sig_gene = NA, err = err_txt2))
    }
    
    # use match() instead of a full logical vector
    ix <- match(gene, rownames(sigdata))
    if (is.na(ix)) {
      .deg_result_png(NULL); .deg_export_df(NULL)
      return(list(deg_sig_gene = NA, err = err_txt2))
    }
    
    row <- sigdata[ix, , drop = FALSE]
    
    # pre-render plot into bytes (tempfile) to avoid device re-runs
    bytes <- .plot_to_raw(function() { print(draw_deg_plot(row, gene)) })
    .deg_result_png(bytes)
    .schedule_ttl(.deg_result_png)
    
    # small export payload only
    .deg_export_df(t(row))
    
    invisible(gc())
    list(deg_sig_gene = row, err = NULL)
  })
  
  output$deg_error2 <- renderUI({
    req(input$main_tab == "DEGs")
    res2 = select_deg()
    if (!is.null(res2$err)) {
      tags$div(style = "color: red; ", res2$err)
    }
  })
  
  output$deg_data_all <- downloadHandler(
    filename = function() {
      paste0(input$cancer,"_",input$datatype,"_PD-L1_DEGs_", Sys.Date(), ".csv")
    },
    content = function(file) {
      res2 <- .get_isolated(deg_cache)   
      df <- res2$deg_sig
      if (is.null(df) || nrow(df) == 0) {
        write.csv(data.frame(Message = "No data to export"), file, row.names = FALSE)
      } else {
        .write_csv_chunked(df, file, chunk_rows = 20000, row_names = TRUE)
        res2$deg_sig <- NULL
        deg_cache(res2)
        invisible(gc())
      }
    }
  )
  
  output$deg_data <- downloadHandler(
    filename = function() {
      paste0(input$cancer,"_",input$datatype,"_PD-L1_DEGs_",input$gene,"_", Sys.Date(), ".csv")
    },
    content = function(file) {
      df <- .get_isolated(.deg_export_df)   
      if (is.null(df)) {
        write.csv(data.frame(Message = "No data to export"), file, row.names = FALSE)
      } else {
        write.csv(df, file, row.names = TRUE)
      }
      .deg_export_df(NULL); invisible(gc())
    }
  )
  
  output$deg <- renderPlot({
    req(input$main_tab == "DEGs")
    res2 <- select_deg()
    data <- res2$deg_sig_gene
    if (all(is.na(data))) return(NULL)
    draw_deg_plot(data, toupper(input$gene))
  })
  
  output$deg_plot <- downloadHandler(
    filename = function() {
      paste0(input$cancer,"_",input$datatype,"_PD-L1_DEGs_",input$gene,"_", Sys.Date(),".png")
    },
    content = function(file) {
      bytes <- .deg_result_png()
      if (is.null(bytes)) {
        res2 <- select_deg()
        if (all(is.na(res2$deg_sig_gene))) return(NULL)
        tf <- tempfile(fileext = ".png")
        on.exit({ if (file.exists(tf)) unlink(tf) }, add = TRUE)
        grDevices::png(tf)
        print(draw_deg_plot(res2$deg_sig_gene, toupper(input$gene)))
        grDevices::dev.off()
        file.copy(tf, file)
      } else {
        writeBin(bytes, file)
      }
      .deg_result_png(NULL); invisible(gc())
    }
  )
  
  
  
  #===========================================================
  #  Correlations
  #===========================================================
  
  observeEvent(list(input$main_tab,input$cancer, input$datatype), {
    req(input$main_tab=="Correlations")
    req(input$cancer,input$datatype)
    
    current_input_cor <- list(a = input$cancer, b = input$datatype)
    test <- is.null(previous_input_cor()) || !identical(current_input_cor, previous_input_cor())
    if (!test) return()
    
    previous_input_cor(current_input_cor)
    cn <- input$cancer; type <- input$datatype
    
    dt3 <- dataprocess1(cn, type)
    expr_data3 <- dt3$exp
    err <- dt3$err
    
    if (is.null(expr_data3)) {
      cor_sig <- NA
      cors <- NA
      cor_cache(list(cors = cors, cor_sig = cor_sig, err = err))
      return(invisible(NULL))
    }
    
    if (type=='Transcriptomic' && (cn %in% rna_noNAT)) {
      expr_t3 <- expr_data3
    } else {
      expr_t3 <- expr_data3[, -grep("_N", colnames(expr_data3)), drop = FALSE]
    }
    
    cors <- cor_dataprocess("CD274", expr_t3)             
    cor_sig <- cors[which(cors$cortype != 'NoSig'), , drop = FALSE]
    cor_cache(list(cors = cors, cor_sig = cor_sig, err = NULL))
    
    .clear_big("dt3", "expr_data3", "expr_t3")
  }, ignoreInit = TRUE)
  
  output$cor_error <- renderUI({
    req(input$main_tab == "Correlations")
    res3 <- cor_cache()
    if (!is.null(res3$err)) {
      tags$div(style = "color: red; font-weight: bold; ", res3$err)
    }
  })
  
  
  .get_cor_sig_now <- function() {
    cn   <- isolate(input$cancer)
    type <- isolate(input$datatype)
    if (is.null(cn) || is.null(type)) return(NULL)
    
    dt  <- dataprocess1(cn, type)
    exp <- dt$exp
    if (is.null(exp)) return(NULL)
    
    if (type == "Transcriptomic" && (cn %in% rna_noNAT)) {
      expr_t <- exp
    } else {
      expr_t <- exp[, -grep("_N", colnames(exp)), drop = FALSE]
    }
    
    cors <- cor_dataprocess("CD274", expr_t)              
    sig  <- cors[which(cors$cortype != "NoSig"), , drop = FALSE]
    rm(dt, exp, expr_t, cors); invisible(gc())
    sig
  }
  
  output$cor_all <- downloadHandler(
    filename = function() {
      paste0(input$cancer, "_", input$datatype, "_PD-L1_Correlations", Sys.Date(), ".csv")
    },
    content = function(file) {
      res3 <- isolate(cor_cache())
      df <- NULL
      if (!is.null(res3) && !is.null(res3$cor_sig) && nrow(res3$cor_sig) > 0) {
        df <- res3$cor_sig
      } else {
        df <- .get_cor_sig_now()
      }
      
      if (is.null(df) || nrow(df) == 0) {
        write.csv(data.frame(Message = "No data to export"), file, row.names = FALSE)
      } else {
        if (exists(".write_csv_chunked", mode = "function")) {
          .write_csv_chunked(df, file, chunk_rows = 20000, row_names = TRUE)
        } else {
          write.csv(df, file, row.names = TRUE)
        }
      }
      invisible(gc())
    }
  )
  
  
  
  ##### select genes
  select_cor <- eventReactive(input$cor_go, {
    gene2 <- toupper(input$gene2)
    
    # 1) read cache without binding dependencies (no re-compute)
    ct <- .get_isolated(cor_cache)
    cordata <- ct$cors
    
    if (is.null(cordata) || is.null(rownames(cordata))) {
      .cor_result_png(NULL); .cor_export_df(NULL)
      return(list(cor_gene = NA, err = err_txt2))
    }
    
    # 2) tiny lookup without building a huge logical vector
    geneset <- c(gene2, "CD274")
    rn <- rownames(cordata)
    ix <- match(geneset, rn)
    ix <- ix[!is.na(ix)]
    if (length(ix) == 0L) {
      .cor_result_png(NULL); .cor_export_df(NULL)
      return(list(cor_gene = NA, err = err_txt2))
    }
    
    cor_gene <- cordata[ix, , drop = FALSE]
    
    # 3) pre-render plot safely into bytes
    bytes <- .plot_to_raw(function() { print(draw_cor_plot(cor_gene, gene2)) })
    .cor_result_png(bytes)
    
    # 4) small export df only
    .cor_export_df(t(cor_gene))
    
    invisible(gc())
    list(cor_gene = cor_gene, err = NULL)
  })
  
  output$cor_error2 <- renderUI({
    req(input$main_tab == "Correlations")
    res3 <- select_cor()
    if (!is.null(res3$err)) {
      tags$div(style = "color: red; ", res3$err)
    }
  })
  
  output$cor_data <- downloadHandler(
    filename = function() {
      paste0(input$cancer,"_",input$datatype,"_PD-L1_Cors_",input$gene,"_", Sys.Date(), ".csv")
    },
    content = function(file) {
      df <- .get_isolated(.cor_export_df)
      if (is.null(df)) {
        write.csv(data.frame(Message = "No data to export"), file, row.names = FALSE)
      } else {
        write.csv(df, file, row.names = TRUE)
      }
      .cor_export_df(NULL); invisible(gc())
    }
  )
  
  output$cor <- renderPlot({
    req(input$main_tab == "Correlations")
    res3 <- select_cor()
    data <- res3$cor_gene
    if (all(is.na(data))) return(NULL)
    draw_cor_plot(data, toupper(input$gene2))
  })
  
  output$cor_plot <- downloadHandler(
    filename = function() {
      paste0(input$cancer,"_",input$datatype,"_PD-L1_Cors_",input$gene,"_", Sys.Date(),".png")
    },
    content = function(file) {
      bytes <- .cor_result_png()
      if (is.null(bytes)) {
        res3 <- select_cor()
        if (all(is.na(res3$cor_gene))) return(NULL)
        tf <- tempfile(fileext = ".png")
        on.exit({ if (file.exists(tf)) unlink(tf) }, add = TRUE)
        grDevices::png(tf)
        print(draw_cor_plot(res3$cor_gene, toupper(input$gene2)))
        grDevices::dev.off()
        file.copy(tf, file)
      } else {
        writeBin(bytes, file)
      }
      .cor_result_png(NULL); invisible(gc())
    }
  )
  
  
  
  #===========================================================
  #  Pathway Activity （ssGSEA）
  #===========================================================
  

  
  observeEvent(list(input$main_tab,input$cancer, input$datatype), {
    req(input$main_tab=="Pathway Activity")
    req(input$cancer,input$datatype)
    
    current_input <- list(a = input$cancer, b = input$datatype)
    test1 <- is.null(previous_input()) || !identical(current_input, previous_input())
    if (!test1) return()
    
    previous_input(current_input)
    
    result <- ssgsea_testprocess(input$cancer, input$datatype)  
    ssgsea_cache(result)
    .ssgsea_t_cache(result$ssgsea_t) 
    invisible(gc())
  }, ignoreInit = TRUE)
  

  path_test <- reactive({
    req(ssgsea_cache())
    res4 <- ssgsea_cache()
    ssgsea_data <- res4$ssgsea_data

    vec <- ssgsea_data[rownames(ssgsea_data)=='PD-L1 expression and PD-1 checkpoint pathway in cancer', , drop = TRUE]
    ssgsea_PD_L1_path <- data.frame(ssgsea_PD_L1_path = as.numeric(vec))
    rownames(ssgsea_PD_L1_path) <- names(vec)
    
    .pdl1_vec(vec)
    .path_export_df(ssgsea_PD_L1_path)
     
    has_normal <- any(grepl("_N", rownames(ssgsea_PD_L1_path)))
    bytes <- .plot_to_raw(function() {
      if (has_normal) {
        print(draw_path_plot(ssgsea_PD_L1_path))
      } else {
        print(draw_path_plot_onlyT(ssgsea_PD_L1_path))
      }
    })
    .path_result_png(bytes); .schedule_ttl(.path_result_png)
    
    invisible(gc())
    ssgsea_PD_L1_path
  })
  
  output$path_error <- renderUI({
    res4 <- ssgsea_cache()
    if (!is.null(res4$err)) {
      tags$div(style = "color: red; ", res4$err)
    }
  })
  
  output$ssgsea_data <- downloadHandler(
    filename = function() {
      paste0(input$cancer,"_",input$datatype,"_PD-L1_pathway_activity_", Sys.Date(), ".csv")
    },
    content = function(file) {
      res4_2 <- .path_export_df()
      write.csv(res4_2, file, row.names = TRUE)
      .path_export_df(NULL); invisible(gc())
    }
  )
  
  output$ssgsea <- renderPlot({
    req(input$main_tab == "Pathway Activity")
    req(input$cancer, input$datatype)
    res4_2 <- path_test()
    has_normal <- any(grepl("_N", rownames(res4_2)))
    if (has_normal) {
      draw_path_plot(res4_2)
    } else {
      draw_path_plot_onlyT(res4_2)
    }
  })
  
  output$ssgsea_plot <- downloadHandler(
    filename = function() {
      paste0(input$cancer,"_",input$datatype,"_PD-L1_pathway_activity_", Sys.Date(),".png")
    },
    content = function(file) {
      bytes <- .path_result_png()
      if (is.null(bytes)) {
        req(input$cancer, input$datatype)
        res4_2 <- path_test()
        tf <- tempfile(fileext = ".png")
        on.exit({ if (file.exists(tf)) unlink(tf) }, add = TRUE)
        grDevices::png(tf)
        has_normal <- any(grepl("_N", rownames(res4_2)))
        if (has_normal) {
          print(draw_path_plot(res4_2))
        } else {
          print(draw_path_plot_onlyT(res4_2))
        }
        grDevices::dev.off()
        file.copy(tf, file)
      } else {
        writeBin(bytes, file)
      }
      .path_result_png(NULL); invisible(gc())
    }
  )
  
  
  
  #===========================================================
  #  Pathway crosstalk
  #===========================================================
  
  output$cs_notice <- renderUI({
    req(input$main_tab == "Pathway crosstalk")
    ssgsea <- ssgsea_cache()
    type <- input$datatype
    if (is.null(ssgsea)) {
      showModal(modalDialog(
        title = "Notice",
        HTML("<span style=' font-size:16px; font-weight:bold;'> 
         Please run the ssGSEA analysis on the [Pathway Activity] page first!
       </span>"),
        easyClose = TRUE
      ))
    } else {
      txt <- paste0('The current cancer type: ', ssgsea$cn,"<br/>",
                    "The current dataset type: ", ssgsea$type,"<br/>",
                    'If you’ve changed the tumor type or dataset, please return to [Pathway Activity] and rerun the analysis.')
      tags$div(style = "color: #4682B4; font-weight: bold;  ", HTML(txt))
    }
  })
  
  path_cs_test <- reactive({
    req(input$main_tab == "Pathway crosstalk")
    req(input$cancer, input$datatype)
    if (!is.null(ssgsea_cache())) {
      st <- .ssgsea_t_cache()
      if (is.null(st)) st <- ssgsea_cache()$ssgsea_t
      if (is.null(st)) return(NULL)                
      ssgsea_t <- data.frame(st, check.names = FALSE)
      
      item <- 'PD-L1 expression and PD-1 checkpoint pathway in cancer'
      ssgsea_cor <- cor_dataprocess(item, ssgsea_t)
      ssgsea_cor_sig <- ssgsea_cor[which(ssgsea_cor$cortype!= 'NoSig'), , drop = FALSE]
      .path_cs_all_df(ssgsea_cor_sig)   # cache for download
      list(ssgsea_cor = ssgsea_cor, ssgsea_cor_sig = ssgsea_cor_sig)
    }
  })
  
  output$path_cs_all <- downloadHandler(
    filename = function() {
      paste0(input$cancer,"_",input$datatype,
             "_PD-L1_pathway_crosstalk_signifcance_", Sys.Date(), ".csv")
    },
    content = function(file) {
      df <- .get_isolated(.path_cs_all_df)

      if (is.null(df) || nrow(df) == 0) {
        res <- isolate(path_cs_test())
        if (!is.null(res)) df <- res$ssgsea_cor_sig
      }
      
      if (is.null(df) || nrow(df) == 0) {
        write.csv(data.frame(Message = "No data to export"), file, row.names = FALSE)
      } else {
        if (exists(".write_csv_chunked", mode = "function")) {
          .write_csv_chunked(df, file, chunk_rows = 20000, row_names = TRUE)
        } else {
          write.csv(df, file, row.names = TRUE)
        }
      }
      invisible(gc())
    }
  )
  
  
  ##### select paths
  select_path_cs_test <- reactive({
    req(input$path_cs)
    res5 <- path_cs_test()
    pathcor_data <- res5$ssgsea_cor
    pathset <- c(input$path_cs, 'PD-L1 expression and PD-1 checkpoint pathway in cancer')
    pathset_cor <- pathcor_data[rownames(pathcor_data) %in% pathset, , drop = FALSE]
    
    bytes <- .plot_to_raw(function() { print(draw_cor_plot(pathset_cor, input$path_cs)) })
    .cs_result_png(bytes); .schedule_ttl(.cs_result_png)

    .cs_export_df(t(pathset_cor))
    pathset_cor
  })
  
  output$cs_data <- downloadHandler(
    filename = function() {
      paste0(input$cancer,"_",input$datatype,"_PD-L1_",input$path_cs,"_crosstalk_", Sys.Date(), ".csv")
    },
    content = function(file) {
      df <- .cs_export_df()
      write.csv(df, file, row.names = TRUE)
      .cs_export_df(NULL); invisible(gc())
    }
  )
  
  output$cs <- renderPlot({
    req(input$path_cs)
    res5 <- select_path_cs_test()
    draw_cor_plot(res5, input$path_cs)
  })
  
  output$cs_plot <- downloadHandler(
    filename = function() {
      paste0(input$cancer,"_",input$datatype,"_PD-L1_",input$path_cs,"_crosstalk_", Sys.Date(),".png")
    },
    content = function(file) {
      bytes <- .cs_result_png()
      if (is.null(bytes)) {
        res5 <- select_path_cs_test()
        tf <- tempfile(fileext = ".png")
        on.exit({ if (file.exists(tf)) unlink(tf) }, add = TRUE)
        grDevices::png(tf)
        print(draw_cor_plot(res5, input$path_cs))
        grDevices::dev.off()
        file.copy(tf, file)
      } else {
        writeBin(bytes, file)
      }
      .cs_result_png(NULL); invisible(gc())
    }
  )
  
  
  #===========================================================
  #  Phosphorylation
  #===========================================================
  
  output$pho_notice <- renderUI({
    req(input$main_tab == "Phosphorylation")
    ssgsea <- ssgsea_cache()
    type <- input$datatype
    if (is.null(ssgsea)) {
      showModal(modalDialog(
        title = "Notice",
        HTML("<span style=' font-size:16px; font-weight:bold;'>
         Please run the ssGSEA analysis on the [Pathway Activity] page first!
         Please use Proteomics data.
       </span>"),
        easyClose = TRUE
      ))
    } else {
      if (type  == "Transcriptomic") {
        showModal(modalDialog(
          title = "Notice",
          HTML("<span style=' font-size:16px; font-weight:bold;'>
         Please rerun the ssGSEA analysis using Proteomics data on the [Pathway Activity] page
       </span>"),
          easyClose = TRUE
        ))
      }
      txt <- paste0('The current cancer type: ', ssgsea$cn,"<br/>",
                    "The current dataset type: ", ssgsea$type,"<br/>",
                    'Please use Proteomics data.', "<br/>",
                    'If you’ve changed the tumor type or dataset, please return to [Pathway Activity] and rerun the analysis.')
      tags$div(style = "color: #4682B4; font-weight: bold;  ",HTML(txt))
    }
  })
  
  observeEvent(input$main_tab, {
    if (identical(input$main_tab, "Phosphorylation")) {
      res <- ssgsea_cache()
      if (!is.null(res) && !is.null(res$ssgsea_data)) {
        pn <- "PD-L1 expression and PD-1 checkpoint pathway in cancer"
        v  <- res$ssgsea_data[pn, , drop = TRUE]
        if (!is.null(colnames(res$ssgsea_data))) names(v) <- colnames(res$ssgsea_data)
        .pdl1_vec(v)
        .ssgsea_t_cache(res$ssgsea_t)  
        res$ssgsea_data <- NULL; res$ssgsea_t <- NULL
        ssgsea_cache(res)
        invisible(gc())
      }
    }
  }, ignoreInit = TRUE)
  
  .ensure_pho_full <- function() {
    df_full <- .pho_full_df()
    if (!is.null(df_full)) return(df_full)
    
    v <- .pdl1_vec(); if (is.null(v)) return(NULL)
    df <- data.frame(ssgsea_PD_L1_path = as.numeric(v))
    rownames(df) <- names(v)
    
    res <- pho_dataprocess(isolate(input$cancer), df) 
    if (!is.null(res) && nrow(res)) {
      .pho_full_df(res)
      .pho_sig_all_df(res[res$change != "NoSignifi", , drop = FALSE])
    }
    invisible(gc())
    .pho_full_df()
  }
  
  pho_test <- reactive({
    req(input$main_tab == "Phosphorylation")
    v <- .pdl1_vec(); if (is.null(v)) return(NULL)
    df <- data.frame(ssgsea_PD_L1_path = as.numeric(v))
    rownames(df) <- names(v)
    cn <- input$cancer
    pho_dataprocess(cn, df)  
  })
  
  
  #### select phosites
  output$pho_inputnotice <- renderUI({
    pho_ex = "Input Example: STAT1_S727"
    tags$div(style = "color: #808080; font-style: italic; font-size:12px;", pho_ex)
  })
  
  
  select_pho <- eventReactive(input$pho_go, {
    site <- toupper(isolate(input$phosite)) 
    res  <- .ensure_pho_full()
    if (is.null(res) || !(site %in% rownames(res))) {
      .pho_result_png(NULL); .pho_export_df(NULL); .pho_site_last(NULL)
      return(NA)
    }
    rowp <- res[site, , drop = FALSE]
    .pho_site_last(site)
    
    bytes <- .plot_to_raw(function() { print(draw_pho_plot(rowp, site)) })
    .pho_result_png(bytes)
    .pho_export_df(t(rowp))  
    
    invisible(gc())
    rowp
  })
  

  select_pho_test <- reactive({
    req(input$main_tab == "Phosphorylation")
    req(input$pho_go)                  
    dt <- select_pho()              

    if (is.null(dt) || (is.atomic(dt) && length(dt) == 1L && is.na(dt))) {
      list(pho = NA, err = err_txt3)
    } else {
      list(pho = dt, err = NULL)
    }
  })
  
  output$pho_error <- renderUI({
    req(input$main_tab == "Phosphorylation")
    res6 <- select_pho_test() 
    if (!is.null(res6$err)) {
      tags$div(style = "color: red; ", res6$err)
    }
  })
  
  output$pho_data <- downloadHandler(
    filename = function() {
      paste0(input$cancer, "_", input$datatype,
             "_PD-L1_High vs Low_", input$phosite, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      df <- .pho_export_df()
      if (is.null(df)) {
        write.csv(data.frame(Message = "No data to export"), file, row.names = FALSE)
      } else {
        write.csv(df, file, row.names = TRUE)
      }
      .pho_export_df(NULL); invisible(gc())
    }
  )
  
  output$pho <- renderPlot({
    req(input$main_tab == "Phosphorylation")
    rowdf <- select_pho()             
    if (all(is.na(rowdf))) return(NULL)
    site <- .pho_site_last(); req(site) 
    draw_pho_plot(rowdf, site)
  })
  
  output$pho_plot <- downloadHandler(
    filename = function() {
      paste0(input$cancer, "_", input$datatype,
             "_PD-L1_High vs Low_", input$phosite, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      bytes <- .pho_result_png()
      if (is.null(bytes)) {
        site <- .pho_site_last(); req(site)
        rowdf <- select_pho(); if (all(is.na(rowdf))) return(NULL)
        tf <- tempfile(fileext = ".png"); on.exit({ if (file.exists(tf)) unlink(tf) }, add = TRUE)
        grDevices::png(tf); print(draw_pho_plot(rowdf, site)); grDevices::dev.off()
        file.copy(tf, file)
      } else {
        writeBin(bytes, file)
      }
      .pho_result_png(NULL); invisible(gc())
    }
  )
  
  output$pho_sig_all <- downloadHandler(
    filename = function() {
      paste0(input$cancer, "_", input$datatype,
             "_PD-L1_pathway_Phosphosites_signifcance_", Sys.Date(), ".csv")
    },
    content = function(file) {
      id <- showNotification("Analyzing phosphosites…", type = "message", duration = NULL)
      
      df <- .pho_sig_all_df()
      if (is.null(df) || nrow(df) == 0L) {
        full <- .ensure_pho_full() 
        if (!is.null(full) && nrow(full)) {
          df <- full[full$change != "NoSignifi", , drop = FALSE]
          .pho_sig_all_df(df)     
        }
      }
      
      if (is.null(df) || nrow(df) == 0L) {
        write.csv(data.frame(Message = "No data to export"), file, row.names = FALSE)
      } else if (exists(".write_csv_chunked", mode = "function")) {
        .write_csv_chunked(df, file, chunk_rows = 20000, row_names = TRUE)
      } else {
        write.csv(df, file, row.names = TRUE)
      }
      
      removeNotification(id)
      invisible(gc())
    }
  )
  
  
}




########################################################################################################
# Final line of code     Final line of code     Final line of code     Final line of code     
########################################################################################################

shinyApp(ui = ui, server = server)