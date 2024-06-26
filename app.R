#### LOAD PACKAGES & FUNCTIONS ####

## First specify the packages of interest

packages <- c("shiny",
              "shinydashboard",
              "tidyr",
              "ggplot2",
              "ggrepel",
              "DT",
              "shinyjs",
              "stringr",
              "RColorBrewer")

## Now load or install&load all
package.check <- lapply(
  
  packages,
  
  FUN = function(x) {
    
    if (!require(x, character.only = TRUE)) {
      
      install.packages(x, dependencies = TRUE)
      
      library(x, character.only = TRUE)
      
    }
    
  }
  
)

# Here define packages which need to be loaded through biocmanager

biocmanager_packages <- c("decoupleR",
                          "DESeq2") 

bioc_package.check <- lapply(
  
  biocmanager_packages,
  
  FUN = function(x) {
    
    if (!require(x, character.only = TRUE)) {
      
      if (!requireNamespace("BiocManager", quietly = TRUE)){
        
        install.packages("BiocManager")
        
      }
      
      BiocManager::install(x, dependencies = TRUE)
      
      library(x, character.only = TRUE)
      
    }
  }
)

# Here define packages to be loaded through Github

github_packages <- c("RAPToR",
                     "wormRef") 

if (!requireNamespace("devtools", quietly = TRUE)){
  install.packages("devtools")}

gh_package.check <- lapply(
  github_packages,
  FUN = function(x) {
    if (!require(str_remove(x, ".*\\/"), character.only = TRUE)) {
      
      if (!requireNamespace("devtools", quietly = TRUE)){
        install.packages("devtools")}
      
      devtools::install_github(x, build_vignettes = TRUE)
      
      library(str_remove(x, ".*\\/"), character.only = TRUE)
      
    }
  }
)

# Define UI for the dashboard
ui <- dashboardPage(
  dashboardHeader(title = "CelEsT: TF activity estimation in C. elegans"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "home", icon = icon("home")),
      menuItem("TF act from DE stats", tabName = "upload_compute", icon = icon("chart-bar")),
      menuItem("DE & TF act from counts", tabName = "alternative_analysis", icon = icon("magnifying-glass-chart"))
    )
  ),
  dashboardBody(
    useShinyjs(),  # Initialize shinyjs
    tabItems(
      # Home tab content
      tabItem(tabName = "home",
              fluidRow(
                column(12, align = "center",
                       img(src = "/logo.png", height = "300px"),
                       h2(HTML("Welcome to <em>Cel</em>EsT")),
                       p(HTML("This tool helps you estimate transcription factor (TF) activity in <em>C. elegans</em>.<br><br>
                        In the tabs on the left, you will find two different analyses available to you:<br><br>
                        
                        1) In 'TF act from DE stats' you can input statistics from an existing differential expression analysis to estimate differential TF activity. <br><br>
                        
                        2) In 'DE & TF act from counts' you can input raw (i.e. not normalised) gene-level RNA-seq read counts. The app will perform a differential expression analysis using <em>DESeq2</em> and then estimate TF activities based on the DE statistics.
                        You can also opt to exclude noise from any effects of your experimental intervention on developmental speed by estimating sample ages using the <em>RAPToR</em> package and incorporating this information in the DE analysis. <br><br>
                        
                        Make sure your data is formatted correctly to ensure accurate results. <em>Cel</em>EsT accepts tab- or space-delimited text files. Gene IDs should be in the first column and ideally be WormBase Gene IDs (e.g. WBGene00020142), sequence names (e.g. T01C8.1) or entrezGene IDs (e.g. 181727). Below you can download sample data for each analysis:<br><br>
                              
                              <a href='sample_data1.txt' download>from DE stats</a> <br>
                              <a href='sample_data2.txt' download>from counts</a> <br><br>
                              
                        This app estimates TF activities using the <em>Cel</em>EsT gene regulatory network with the <em>decoupleR</em> package, employing the multivariate linear model ('mlm') method.<br><br>

                        If you use <em>Cel</em>EsT in a publication, <b>don't forget to cite:</b><br><br>
                        
                        Perez, MF (2024) CelEsT: a unified gene regulatory network for estimating transcription factor activity in C. elegans. <em>bioRxiv</em> <a href='https://doi.org/10.1093/bioadv/vbac016' download>doi: 10.1093/bioadv/vbac016</a> <br><br> & <br><br>
                        
                        Badia-i-Mompel, P et al. (2022) decoupleR: ensemble of computational methods to infer biological activities from omics data. <em>Bioinformatics Advances</em>, <a href='https://doi.org/10.1093/bioadv/vbac016' download>doi: 10.1093/bioadv/vbac016</a> <br>     <br><br>
                        
                        If you use this app for DE analysis, <b>cite:</b><br><br>
                        
                        Love, MI et al. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. <em>Genome Biology</em>, <a href='https://doi.org/10.1186/s13059-014-0550-8' download>doi: 10.1186/s13059-014-0550-8</a><br><br>
                        
                        If you choose to control for developmental age in your analysis, <b>cite:</b><br><br>
                        
                        Bulteau, R & Francesconi, M. (2022) Real age prediction from the transcriptome with RAPToR. <em>Nature Methods</em>, <a href='https://doi.org/10.1038/s41592-022-01450-0' download>doi: 10.1038/s41592-022-01450-0</a>
                        
                              "))
                )
              )
      ),
      # Upload & Compute tab content
      tabItem(tabName = "upload_compute",
              fluidRow(
                shinydashboard::box(
                  title = "TF activity from DE stats",
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  width = 12,
                  fileInput("userdata", "Input DE_stats", placeholder = "Two-column tab-separated file with gene IDs in 1st column and DE stats in 2nd column"),
                  actionButton("compute", "Compute TF activities", class = "btn-primary"),
                  textOutput("status_message")
                )
              ),
              actionButton("reset_button", "Reset", style = "display: none;"),
              fluidRow(
                uiOutput("volcano_ui")
              ),
              fluidRow(
                uiOutput("table_ui")
              ),
              div(id = "tab1status", style = "text-align: center; display: none;",
                  tags$h4("The analysis will take a few minutes. Do not navigate away from this page.Current step: Analysis initiated"))
      ),
      tabItem(tabName = "alternative_analysis",
              fluidRow(
                shinydashboard::box(
                  title = "DE analysis & TF activity from RNA-seq raw counts",
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  width = 12,
                  fileInput("counts_userdata", "Input Data", placeholder = "Tab-separated file with gene IDs in first column and sample names as column headers"),
                  tableOutput("column_names"),
                  textInput("control_samples", "Control Samples", placeholder = "Enter control sample names separated by commas (without quotation marks)"),
                  textInput("treatment_samples", "Treatment Samples", placeholder = "Enter treatment sample names separated by comma (without quotation marks)s"),
                  checkboxInput("correct_speed", "Correct for developmental speed using RAPToR?", value = FALSE),
                  uiOutput("age_selection"),
                  actionButton("counts_compute", "Compute DE & TF activities", class = "btn-primary")
                )
              ),
              actionButton("reset_button2", "Reset", style = "display: none;"),
              fluidRow(
                uiOutput("counts_volcano_ui")
              ),
              fluidRow(
                uiOutput("counts_table_ui")
              ),
              fluidRow(
                uiOutput("counts_DE_ui")
              ),
              fluidRow(
                uiOutput("age_DT_box")
              ),
              div(id = "tab2status", style = "text-align: center; display: none;",
                  tags$h4("The analysis will take a few minutes. Do not navigate away from this page.")
              )
      )
    )
  )
)

server <- function(input, output, session) {
  
  counts_error_occurred <- reactiveVal(FALSE)
  fromDE_error_occurred <- reactiveVal(FALSE)
  
  control_samples <- reactiveVal(NULL)
  treatment_samples <- reactiveVal(NULL)
  
  map2color <- function(x, 
                        pal,
                        limits = NULL){
    
    if(is.null(limits)){
      limits = range(x)
    }
    
    pal[findInterval(x, seq(limits[1], limits[2], length.out = length(pal) + 1), all.inside = TRUE)]
    
  }
  
  CelEsT <- read.table("./www/CelEsT_GRN.txt",
                       sep = "\t",
                       header = TRUE)
  
  CelEsT_BM <- readRDS("./www/CelEsT_BM.rds")
  
  fullset_TFs_BM <- read.table("./www/fullset_TFs_BM_with_MOR_for_TS2.txt",
                               header = TRUE,
                               sep = "\t")

#### CODE FOR ANALYSIS FROM DE STATS ####
  
  
  # Reactive expression to handle data upload
  uploading <- reactive({
    req(input$userdata)
    userdata <- read.table(input$userdata$datapath, header = TRUE)
    userdata[is.na(userdata)] <- 0
    userdata
  })
  
  clearText <- observeEvent(input$compute, {
    output$status_message <- renderText(NULL)
  })
  
  observeEvent(input$reset_button, {
    # Reload the session
    session$reload()
    runjs('document.getElementById("tab1status").style.display = "none";')
    
  })
  
  observeEvent(input$reset_button2, {
    # Reload the session
    session$reload()
    runjs('document.getElementById("tab2status").style.display = "none";')
    
  })
  
  # Reactive expression to handle computation
  computation <- eventReactive(input$compute, {
    
    runjs('document.getElementById("tab1status").style.display = "block";')

    tryCatch({
    
    userdata <- uploading()
    
    runjs(paste0('document.getElementById("tab1status").innerHTML = "', "The analysis will take a few minutes. Do not navigate away from this page. Current step: <b>TF computation has begun</b>", '";'))

    shiny::validate(need(ncol(userdata) == 2, "Data must be 2 columns; the first for gene IDs and one for DE stats"))
    
    # restrict to genes that we have in CelEsT
    userdata <- userdata[str_remove(userdata[, 1], "\\.[a-z]{1,2}") %in% unlist(CelEsT_BM), ]
    
    # setting gene names directly will give error in case of duplications.. 
    geneids_fromuserdata <- userdata[, 1]
    
    if(any(duplicated(geneids_fromuserdata))){
      
      nmbr_dupl <- sum((duplicated(geneids_fromuserdata)|duplicated(geneids_fromuserdata, fromLast = TRUE)))
      
      uniquegeneids_fromuserdata <- geneids_fromuserdata[!(duplicated(geneids_fromuserdata)|duplicated(geneids_fromuserdata, fromLast = TRUE))]
      userdata <- userdata[match(uniquegeneids_fromuserdata, userdata$geneid), ]
      
      row.names(userdata) <- uniquegeneids_fromuserdata
      
      showModal(modalDialog(
        title = "Duplicated Genes Removed",
        paste0(nmbr_dupl, " genes have been removed from analysis due to duplicated gene IDs"),
        easyClose = TRUE
      ))

      
    } else {
      
      row.names(userdata) <- geneids_fromuserdata
      
    }
    
    userdata <- userdata[, 2, drop = FALSE]
    
    userdata[is.na(userdata)] <- 0
    
    gseq_check <- any(str_remove(row.names(userdata), "\\.[a-z]{1,2}$") %in% CelEsT_BM$wormbase_gseq)
    entrezGeneID_check <- any(row.names(userdata) %in% CelEsT_BM$entrezgene_id)
    WBgeneID_check <- any(row.names(userdata) %in% CelEsT_BM$wormbase_gene)
    genename_check <- any(row.names(userdata) %in% CelEsT_BM$external_gene_name)
    
    shiny::validate(need(sum(c(gseq_check, entrezGeneID_check, WBgeneID_check, genename_check)) == 1,
                  "Gene IDs are inconsistent or absent; please revise. We accept WormBase sequence IDs, entrezgene IDs or WormBase IDs"))
    
    if(isFALSE(gseq_check)){
      
      if(isTRUE(entrezGeneID_check)){
        
        convertfromentrez <- CelEsT_BM[match(row.names(userdata), CelEsT_BM$entrezgene_id), c("wormbase_gseq", "entrezgene_id")]
        uniquegseq_fromentrez <- convertfromentrez[!(duplicated(convertfromentrez$wormbase_gseq)|duplicated(convertfromentrez$wormbase_gseq, fromLast = TRUE)), "wormbase_gseq"]
        userdata <- userdata[match(convertfromentrez[match(uniquegseq_fromentrez, convertfromentrez$wormbase_gseq), "entrezgene_id"], row.names(userdata)), , drop = FALSE]
        
        row.names(userdata) <- uniquegseq_fromentrez
        
      }
      
      if(isTRUE(WBgeneID_check)){
        
        convertfromWBID <- CelEsT_BM[match(row.names(userdata), CelEsT_BM$wormbase_gene), c("wormbase_gseq", "wormbase_gene")]
        uniquegseq_fromWBID <- convertfromWBID[!(duplicated(convertfromWBID$wormbase_gseq)|duplicated(convertfromWBID$wormbase_gseq, fromLast = TRUE)), "wormbase_gseq"]
        userdata <- userdata[match(convertfromWBID[match(uniquegseq_fromWBID, convertfromWBID$wormbase_gseq), "wormbase_gene"], row.names(userdata)), , drop = FALSE]
        
        row.names(userdata) <- uniquegseq_fromWBID
        
      }
      
    } else {
      
      # remove trailing letter from CDS id if applicable
      row.names(userdata) <- str_remove(row.names(userdata), "\\.[a-z]{1,2}$")
      
    }
    
    userdata_decouple <- decoupleR::decouple(
      mat = userdata[, 1, drop = FALSE], 
      network = CelEsT,
      .source = "source",
      .target = "target",
      statistics = "mlm",
      args = list(mlm = list(.mor = "weight")),
      consensus_score = FALSE
    )
    
    # Hide the notification
    runjs('document.getElementById("tab1status").style.display = "none";')
    
    userdata_decouple
    
    }, error = function(e){
      
      output$status_message <- renderText({
        
        paste("An error occurred:", e$message)

      })
      
      fromDE_error_occurred(TRUE)
      # Hide the notification
      runjs('document.getElementById("tab1status").style.display = "none";')
      return(NULL)
      
    })
    
  })
  
  # Reactive expression to prepare data for plot
  userdata_forplot <- reactive({
    
    userdata_decouple <- computation()
    
    userdata_decouple_wide <- pivot_wider(
      data = userdata_decouple[, c("source", "condition", "score")], 
      names_from = "source", 
      values_from = "score"
    )
    userdata_decouple_wide <- data.frame(t(userdata_decouple_wide))
    
    colnames(userdata_decouple_wide) <- userdata_decouple_wide[1, ]
    userdata_decouple_wide <- userdata_decouple_wide[-1, , drop = FALSE]
    
    userdata_decouple_wide[] <- lapply(userdata_decouple_wide, as.numeric)
    
    userdata_decouple_wide_p <- as.data.frame(
      pivot_wider(data = userdata_decouple[, c("source", "condition", "p_value")], 
                  names_from = "source", values_from = "p_value")
    )
    row.names(userdata_decouple_wide_p) <- unlist(userdata_decouple_wide_p[, 1])
    userdata_decouple_wide_p <- userdata_decouple_wide_p[, -1, drop = FALSE]
    
    userdata_decouple_wide_p <- data.frame(t(userdata_decouple_wide_p))
    
    userdata_decouple_wide_p[] <- lapply(userdata_decouple_wide_p, as.numeric)
    
    userdata_forbubbleplot <- data.frame(
      "mean_z" = unlist(userdata_decouple_wide[row.names(userdata_decouple_wide_p), ]),
      "p_val" = -log10(unlist(userdata_decouple_wide_p)),
      "label" = unlist(CelEsT_BM[match(row.names(userdata_decouple_wide_p), CelEsT_BM$wormbase_gseq), "usefullabel"])
    )
    
    userdata_forbubbleplot[userdata_forbubbleplot$p_val < 2.2, "label"] <- ""
    
    userdata_forbubbleplot
    
  })
  
  userdata_forDT <- reactive({
    userdata_decouple <- computation()
    
    userdata_decouple_wide <- pivot_wider(
      data = userdata_decouple[, c("source", "condition", "score")], 
      names_from = "source", 
      values_from = "score"
    )
    
    userdata_decouple_wide <- data.frame(t(userdata_decouple_wide))
    
    colnames(userdata_decouple_wide) <- userdata_decouple_wide[1, ]
    userdata_decouple_wide <- userdata_decouple_wide[-1, , drop = FALSE]
    
    userdata_decouple_wide[] <- lapply(userdata_decouple_wide, as.numeric)
    
    userdata_decouple_wide_p <- as.data.frame(
      pivot_wider(data = userdata_decouple[, c("source", "condition", "p_value")], 
                  names_from = "source", values_from = "p_value")
    )
    row.names(userdata_decouple_wide_p) <- unlist(userdata_decouple_wide_p[, 1])
    userdata_decouple_wide_p <- userdata_decouple_wide_p[, -1, drop = FALSE]
    
    userdata_decouple_wide_p <- data.frame(t(userdata_decouple_wide_p))
    
    userdata_decouple_wide_p[] <- lapply(userdata_decouple_wide_p, as.numeric)
    
    userdata_forDT <- data.frame(
      "Score" = unlist(userdata_decouple_wide[row.names(userdata_decouple_wide_p), ]),
      "Pval" = unlist(userdata_decouple_wide_p),
      "Sequence_name" = row.names(userdata_decouple_wide_p),
      "Gene_name" = unlist(CelEsT_BM[match(row.names(userdata_decouple_wide_p), CelEsT_BM$wormbase_gseq), "wormbase_locus"]),
      "Entrez_gene_ID" = unlist(CelEsT_BM[match(row.names(userdata_decouple_wide_p), CelEsT_BM$wormbase_gseq), "entrezgene_id"]),
      "WormBase_MOR" = unlist(fullset_TFs_BM[match(row.names(userdata_decouple_wide_p), fullset_TFs_BM$wormbase_gseq), "Wormbase_mode_of_regulation"]),
      "UniProt_MOR" = unlist(fullset_TFs_BM[match(row.names(userdata_decouple_wide_p), fullset_TFs_BM$wormbase_gseq), "Uniprot_mode_of_regulation"]),
      row.names = NULL
      
    )
    
    userdata_forDT <- userdata_forDT[order(userdata_forDT$Pval), ]
    
    userdata_forDT[, c("WormBase_MOR", "UniProt_MOR")] <- str_replace_all(unlist(userdata_forDT[, c("WormBase_MOR", "UniProt_MOR")]), "-1", "REPRESSOR")
    userdata_forDT[, c("WormBase_MOR", "UniProt_MOR")] <- str_replace_all(unlist(userdata_forDT[, c("WormBase_MOR", "UniProt_MOR")]), "0", "BIFUNCTIONAL")
    userdata_forDT[, c("WormBase_MOR", "UniProt_MOR")] <- str_replace_all(unlist(userdata_forDT[, c("WormBase_MOR", "UniProt_MOR")]), "1", "ACTIVATOR")
    
    userdata_forDT
    
  })
  
  output$volcano_ui <- renderUI({
    req(computation())  # Render the UI only if the computation has been done
    shinydashboard::box(
      title = "TF activity volcano plot",
      status = "primary",
      solidHeader = TRUE,
      collapsible = TRUE,
      width = 12,
      plotOutput("TF_volcano", height = "500px")
    )
  })
  
  output$table_ui <- renderUI({
    req(computation())  # Render the UI only if the computation has been done
    shinydashboard::box(
      title = "TF activity data table",
      status = "primary",
      solidHeader = TRUE,
      collapsible = TRUE,
      width = 12,
      downloadButton("downloadData", "Save Data Table", class = "btn-success"),
      DTOutput("data_table")
    )
  })
  
  output$TF_volcano <- renderPlot({
    data <- userdata_forplot()
    #     data <- userdata_forbubbleplot
    
    tempcolours_for_plot <- map2color(data$p_val, pal = colorRampPalette(c("grey", "grey", "grey", "red", "red", "red"))(100))
    
    ggplot(data, aes(x = mean_z, y = p_val, size = 10^p_val, label = toupper(label))) + 
      geom_point(col = tempcolours_for_plot, alpha = 0.7) + 
      theme_classic() + 
      geom_label_repel(size = 2, max.overlaps = 40) +
      geom_vline(xintercept = 0, linetype = "dashed", col = "grey") + 
      theme(
        legend.position = "none",
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)
      ) + 
      ylab(substitute("-log"[10]~"(p-value)")) + 
      xlab("TF activity") 
  })
  
  # Render the DataTable
  output$data_table <- renderDT({
    data <- userdata_forDT()
    datatable(data,
              rownames = FALSE)
  }, options = list(pageLength = 10))  # You can adjust the number of rows displayed per page
  
  # Download handler for the DataTable
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("TFactivities_fromDE-", Sys.Date(), ".txt", sep = "")
    },
    content = function(file) {
      write.table(userdata_forDT(), file, sep = "\t", row.names = FALSE, col.names = TRUE)
    }
  )
  
  observe({
    if (fromDE_error_occurred()) {
      shinyjs::show("reset_button")
    } else {
      hide("reset_button")
    }
  })
  
#### CODE FOR ANALYSIS FROM COUNTS ####

  observeEvent(input$counts_compute, {
    # Use isolate to update the reactiveVal only when the button is clicked
    control_samples <- isolate(input$control_samples)
    treatment_samples <- isolate(input$treatment_samples)
    
  })
  
  observe({
    
    if (input$correct_speed) {
      
      output$age_selection <- renderUI({
        radioButtons("expected_age", "Expected Age",
                     choices = list("Embryo" = "Cel_embryo",
                                    "Larvae" = "Cel_larval",
                                    "Larvae / Young Adult" = "Cel_larv_YA",
                                    "Day 1/2 Adult" = "Cel_YA_2",
                                    "Day 3+ Adult" = "Cel_aging_1"),
                     selected = "embryo")
      })
    } else {
      
      output$age_selection <- renderUI({
        NULL
      })
    }
  })
  
  
  observeEvent(input$counts_compute, {
    if (is.null(input$counts_userdata)) {
      runjs(paste0('document.getElementById("tab2status").innerHTML = "', "Please input a file before pressing compute", '";'))
    } else{
      runjs('document.getElementById("tab2status").style.display = "block";')
    }
  })
  
  uploading_counts <- reactive({
    req(input$counts_userdata)
    
    userdata <- read.table(input$counts_userdata$datapath, header = TRUE)

    userdata[is.na(userdata)] <- 0
    userdata
  })
  
  output$column_names <- renderPrint({
    req(uploading_counts())
    colnames(uploading_counts())[2:ncol(uploading_counts())]
  })
  
  counts_data_process <- eventReactive(input$counts_compute, {
    
    control_samples <- isolate(input$control_samples)
    control_samples <- str_remove_all(control_samples, " ")
    control_samples <- unlist(str_split(control_samples, ","))
    
    treatment_samples <- isolate(input$treatment_samples)
    treatment_samples <- str_remove_all(treatment_samples, " ")
    treatment_samples <- unlist(str_split(treatment_samples, ","))
    
    # setting gene names directly will give error in case of duplications.. 
    uploading_counts <- uploading_counts()

    geneids_fromuserdata <- uploading_counts[, 1]
    
    if(any(duplicated(geneids_fromuserdata))){
      
      nmbr_dupl <- sum((duplicated(geneids_fromuserdata)|duplicated(geneids_fromuserdata, fromLast = TRUE)))
      
      uniquegeneids_fromuserdata <- geneids_fromuserdata[!(duplicated(geneids_fromuserdata)|duplicated(geneids_fromuserdata, fromLast = TRUE))]
      uploading_counts <- uploading_counts[match(uniquegeneids_fromuserdata, uploading_counts$geneid), ]
      
      row.names(uploading_counts) <- uniquegeneids_fromuserdata
      
      showModal(modalDialog(
        title = "Duplicated Genes Removed",
        paste0(nmbr_dupl, " genes have been removed from analysis due to duplicated gene IDs"),
        easyClose = TRUE
      ))
      
    } else {
      
      row.names(uploading_counts) <- geneids_fromuserdata
      
    }
    
    uploading_counts <- uploading_counts[, 2:ncol(uploading_counts), drop = FALSE]
    
    # replace NAs with 0s
    uploading_counts[is.na(uploading_counts)] <- 0
    
    gseq_check <- any(str_remove(row.names(uploading_counts), "\\.[a-z]{1,2}$") %in% CelEsT_BM$wormbase_gseq)
    entrezGeneID_check <- any(row.names(uploading_counts) %in% CelEsT_BM$entrezgene_id)
    WBgeneID_check <- any(row.names(uploading_counts) %in% CelEsT_BM$wbps_gene_id)
    genename_check <- any(row.names(uploading_counts) %in% CelEsT_BM$external_gene_name)
    
    shiny::validate(need(sum(c(gseq_check, entrezGeneID_check, WBgeneID_check, genename_check)) == 1,
                  "Gene IDs are inconsistent or absent; please revise. We accept WormBase sequence IDs, entrezgene IDs or WormBase IDs"))
    
    if(isFALSE(gseq_check)){
      
      if(isTRUE(entrezGeneID_check)){
        
        convertfromentrez <- CelEsT_BM[match(row.names(uploading_counts), CelEsT_BM$entrezgene_id), c("wormbase_gseq", "entrezgene_id")]
        uniquegseq_fromentrez <- convertfromentrez[!(duplicated(convertfromentrez$wormbase_gseq)|duplicated(convertfromentrez$wormbase_gseq, fromLast = TRUE)), "wormbase_gseq"]
        uploading_counts <- uploading_counts[match(convertfromentrez[match(uniquegseq_fromentrez, convertfromentrez$wormbase_gseq), "entrezgene_id"], row.names(uploading_counts)), , drop = FALSE]
        
        row.names(uploading_counts) <- uniquegseq_fromentrez
        
      }
      
      if(isTRUE(WBgeneID_check)){
        
        convertfromWBID <- CelEsT_BM[match(row.names(uploading_counts), CelEsT_BM$wbps_gene_id), c("wormbase_gseq", "wbps_gene_id")]
        uniquegseq_fromWBID <- convertfromWBID[!(duplicated(convertfromWBID$wormbase_gseq)|duplicated(convertfromWBID$wormbase_gseq, fromLast = TRUE)), "wormbase_gseq"]
        uploading_counts <- uploading_counts[match(convertfromWBID[match(uniquegseq_fromWBID, convertfromWBID$wormbase_gseq), "wbps_gene_id"], row.names(uploading_counts)), , drop = FALSE]
        
        row.names(uploading_counts) <- uniquegseq_fromWBID
        
      }
      
    } else {
      
      # remove trailing letter from CDS id if applicable
      row.names(uploading_counts) <- str_remove(row.names(uploading_counts), "\\.[a-z]{1,2}$")
      
    }
    
    uploading_counts
    
  })
  
  tempref <- eventReactive(input$counts_compute, {
    
    if (input$correct_speed) {
      
      req(input$expected_age)

      tempref <- prepare_refdata(input$expected_age, 'wormRef', 5000)

      return(tempref)
      
    } else {
      
      return(NULL)
      
    }})
  
  agecorr_process <- eventReactive(input$counts_compute, {
    if (input$correct_speed) {
      
      req(input$expected_age)
      req(tempref())
      
      counts_to_tpm <- function(counts, featureLength, meanFragmentLength) {
        
        stopifnot(length(featureLength) == nrow(counts))
        stopifnot(length(meanFragmentLength) == ncol(counts))
        
        effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
          featureLength - meanFragmentLength[i] + 1
        }))
        
        idx <- apply(effLen, 1, function(x) min(x) > 1)
        temp_counts <- counts[idx,]
        temp_effLen <- effLen[idx,]
        temp_featureLength <- featureLength[idx]
        
        tpm <- do.call(cbind, lapply(1:ncol(temp_counts), function(i) {
          rate = log(temp_counts[,i]) - log(temp_effLen[,i])
          denom = log(sum(exp(rate)))
          exp(rate - denom + log(1e6))
        }))
        
        colnames(tpm) <- colnames(temp_counts)
        rownames(tpm) <- rownames(temp_counts)
        
        return(tpm)
        
      }
      
      temp_data <- counts_data_process()  

      control_samples <- isolate(input$control_samples)
      control_samples <- str_remove_all(control_samples, " ")
      control_samples <- unlist(str_split(control_samples, ","))
      
      treatment_samples <- isolate(input$treatment_samples)
      treatment_samples <- str_remove_all(treatment_samples, " ")
      treatment_samples <- unlist(str_split(treatment_samples, ","))
      
      # Show the notification in the center
      runjs(paste0('document.getElementById("tab2status").innerHTML = "', "The analysis will take a few minutes. Do not navigate away from this page. Current step: <b>Preparing to estimate sample ages</b>", '";'))
      
      if (any(!colnames(temp_data)[2:ncol(temp_data)] %in% c(control_samples, treatment_samples))) {

        showNotification("Warning: additional columns are present in uploaded data", type = "warning", duration = NULL)
        
      }
      
      my_temp_lengths <- sapply(row.names(temp_data), function(x){
        
        mean(wormRef::Cel_genes[wormRef::Cel_genes$sequence_name == x, "transcript_length"])
        
      })
      
      my_temp_lengths <- my_temp_lengths[!is.na(my_temp_lengths)]
      
      gene_there <- base::intersect(names(my_temp_lengths), row.names(temp_data))
      
      temp_data_TPM <- counts_to_tpm(counts = temp_data[gene_there, ],
                                     featureLength = my_temp_lengths[gene_there],
                                     meanFragmentLength = as.numeric(rep(100, times = ncol(temp_data))))
      
      temp_data_TPM_expressed <- temp_data_TPM[apply(temp_data_TPM, 1, function(x){any(x != 0)}), ]
      
      temp_data_TPM_expressed_norm <- limma::normalizeBetweenArrays(temp_data_TPM_expressed, method = "quantile")
      temp_data_TPM_expressed_norm_log <- log1p(temp_data_TPM_expressed_norm) # log1p(x) = log(x + 1)
      
      row.names(temp_data_TPM_expressed_norm_log) <- wormRef::Cel_genes[match(row.names(temp_data_TPM_expressed_norm_log), wormRef::Cel_genes$sequence_name), "wb_id"]
      
      runjs(paste0('document.getElementById("tab2status").innerHTML = "', "The analysis will take a few minutes. Do not navigate away from this page. Current step: <b>Estimating sample ages</b>", '";'))
      
      temp_sample_ae <- RAPToR::ae(samp = temp_data_TPM_expressed_norm_log,                         
                                   refdata = tempref())
      
      temp_sample_ae$age.estimates[c(control_samples, treatment_samples), ]
      
    } else {
      
      # Return NULL if age correction is not required
      return(NULL)
    }
  })
  
  counts_DE <- reactive({
    req(input$counts_compute)
    
    if (input$correct_speed) {

      age_estimates <- req(agecorr_process())
    
      runjs(paste0('document.getElementById("tab2status").innerHTML = "', "The analysis will take a few minutes. Do not navigate away from this page. Current step: <b>Performing DE analysis </b>", '";'))
      
      tryCatch({
        req(agecorr_process())
        
        req(tempref())
        
        tempref <- tempref()
        
        temp_data <- counts_data_process()
        
        control_samples <- input$control_samples
        control_samples <- str_remove_all(control_samples, " ")
        control_samples <- unlist(str_split(control_samples, ","))
        
        treatment_samples <- input$treatment_samples
        treatment_samples <- str_remove_all(treatment_samples, " ")
        treatment_samples <- unlist(str_split(treatment_samples, ","))
        
        shiny::validate(need(c(control_samples, treatment_samples) %in% colnames(temp_data), "Control/treatment samples do not match column names"))

        temp_data  <- temp_data [, c(control_samples, treatment_samples)]
        
        if(input$expected_age == "Cel_embryo"){
          interpolation_offset <- 60
        } else {
          interpolation_offset <- 1
        }
        
        interpolation_range <- range(age_estimates[, 1]) + c(-interpolation_offset, interpolation_offset)
        
        interpolation_index <- get_refTP(tempref, ae_values = interpolation_range)
        interpolation_index <- interpolation_index[1]:interpolation_index[2]
        
        interpolated_ref_expression <- get_refTP(tempref,
                                                 ae_values = interpolation_range,
                                                 return.idx = FALSE)
        
        interp_time <- tempref$time[interpolation_index]
        interp_tpm <- tempref$interpGE[, interpolation_index]
        
        wormRef_gene_lengths <- sapply(row.names(interp_tpm), function(x){
          
          mean(wormRef::Cel_genes[wormRef::Cel_genes$wb_id == x, "transcript_length"])
          
        })
        
        libsize <- 25e6
        
        interp_count <- t((t (exp(interp_tpm) - 1)/ 1e6) * (libsize / median(wormRef_gene_lengths))) * wormRef_gene_lengths
        interp_count[interp_count < 0] <- 0
        interp_count <- round(interp_count)
        
        df_SSQ <- sapply(1:8, function(df_param){
          
          sum(residuals(lm(t(interp_tpm) ~ splines::ns(interp_time, df = df_param))) ^2)
          
        })

        threshold <- 0.01
        diff_1 <- diff(df_SSQ)
        
        temp_plateau <- which.max((diff_1 / diff_1[1]) < threshold) 
        
        temp_data <- temp_data[apply(temp_data, 1, max) > 5, ]
        
        uniquegeneids_fromtempdata <- row.names(temp_data)
        uniqueIDs_toWBID <- wormRef::Cel_genes[match(uniquegeneids_fromtempdata, wormRef::Cel_genes$sequence_name), c("wb_id", "sequence_name")]
        
        uniqueIDs_toWBID <- uniqueIDs_toWBID[!(duplicated(uniqueIDs_toWBID$wb_id)|duplicated(uniqueIDs_toWBID$wb_id, fromLast = TRUE)), ]
        
        temp_data <- temp_data[row.names(temp_data) %in% uniqueIDs_toWBID$sequence_name, ]
        
        row.names(temp_data) <- wormRef::Cel_genes[match(row.names(temp_data), wormRef::Cel_genes$sequence_name), "wb_id"]
        
        combine_count <- do.call(cbind, format_to_ref(temp_data, interp_count)[1:2])
        
        tempcoldata <- data.frame(group = c(rep("control", times = length(control_samples)), rep("treatment", times = length(treatment_samples))))
        row.names(tempcoldata) <- c(control_samples, treatment_samples)
        
        combine_coldata <- data.frame(time = c(age_estimates[, 1], interp_time),
                                      strain = c(tempcoldata[row.names(age_estimates), 1], rep("control", ncol(interp_count))),
                                      batch = rep(c("sample", "reference"), c(ncol(temp_data), ncol(interp_count))))
        
        combine_coldata$strain <- factor(combine_coldata$strain, levels = c("control", "treatment"))
        combine_coldata$batch <- factor(combine_coldata$batch, levels = c("sample", "reference"))

        row.names(combine_coldata) <- c(row.names(age_estimates), colnames(interp_count))
        
        dd0 <- DESeqDataSetFromMatrix(countData = combine_count[, 1:length(c(control_samples, treatment_samples))],
                                      colData = combine_coldata[1:length(c(control_samples, treatment_samples)), ],
                                      design = ~ strain)
        
        dd0 <- estimateSizeFactors(dd0)
        dd0 <- estimateDispersions(dd0, fitType = "local") 
        d0 <- dispersions(dd0) # store dispersions
        d0[is.na(d0)] <- 0 # remove NAs
        
        formula <- paste0("~ splines::ns(time, df = ", temp_plateau, ") + batch + strain")

        dd1 <-  DESeqDataSetFromMatrix(
          countData = combine_count,
          colData = combine_coldata,
          design = as.formula(formula) # use df found above
        )
        
        dd1 <- estimateSizeFactors(dd1)

        dispersions(dd1) <- d0
        
        dd1 <- nbinomWaldTest(dd1)
        
        result <- as.data.frame(results(dd1, contrast = c("strain", "treatment", "control")))
        
        return(result)
      }, error = function(e){

        runjs(paste0('document.getElementById("tab2status").innerHTML = "', e$message, '";'))
        counts_error_occurred(TRUE)
        
        
      })
    } else {

      runjs(paste0('document.getElementById("tab2status").innerHTML = "', "The analysis will take a few minutes. Do not navigate away from this page. Current step: <b>Performing DE analysis</b>", '";'))
      
      tryCatch({

        temp_data <- counts_data_process()

        control_samples <- isolate(input$control_samples)
        control_samples <- str_remove_all(control_samples, " ")
        control_samples <- unlist(str_split(control_samples, ","))
        
        treatment_samples <- isolate(input$treatment_samples)
        treatment_samples <- str_remove_all(treatment_samples, " ")
        treatment_samples <- unlist(str_split(treatment_samples, ","))
        
        shiny::validate(need(c(control_samples, treatment_samples) %in% colnames(temp_data), "Control/treatment samples do not match column names"))
        
        if (any(!colnames(temp_data)[2:ncol(temp_data)] %in% c(control_samples, treatment_samples))) {

          showNotification("Warning: additional columns are present in uploaded data", type = "warning", duration = NULL)
          
        }
        
        temp_data <- temp_data[, c(control_samples, treatment_samples)]
        
        tempcoldata <- data.frame(group = c(rep("control", times = length(control_samples)), rep("treatment", times = length(treatment_samples))))
        row.names(tempcoldata) <- c(control_samples, treatment_samples)
        
        tempdds <- DESeqDataSetFromMatrix(countData = temp_data, colData = tempcoldata, design = ~ group)
        tempdds <- DESeq(tempdds)
        
        result <- as.data.frame(results(tempdds))
        
        return(result)
      }, error = function(e){
        runjs(paste0('document.getElementById("tab2status").innerHTML = "', e$message, '";'))
        counts_error_occurred(TRUE)
        return("An error occurred")
      })
    }
  })
  

  counts_computation <- reactive({
    
    req(counts_DE())
    
    if(is.character(counts_DE())){

    } else {

    runjs(paste0('document.getElementById("tab2status").innerHTML = "', "The analysis will take a few minutes. Do not navigate away from this page. Current step: <b>TF activity computation started</b>", '";'))
    }
      
    DEdata <- counts_DE()
    
    control_samples <- input$control_samples
    control_samples <- str_remove_all(control_samples, " ")
    control_samples <- unlist(str_split(control_samples, ","))
    
    treatment_samples <- input$treatment_samples
    treatment_samples <- str_remove_all(treatment_samples, " ")
    treatment_samples <- unlist(str_split(treatment_samples, ","))
    
    DEdata <- DEdata[str_remove(row.names(DEdata), "\\.[a-z]{1,2}") %in% unlist(CelEsT_BM), ]
    
    geneids_fromDEdata <- DEdata[, 1]
    
    DEdata[is.na(DEdata)] <- 0
    
    gseq_check <- any(str_remove(row.names(DEdata), "\\.[a-z]{1,2}$") %in% CelEsT_BM$wormbase_gseq)
    entrezGeneID_check <- any(row.names(DEdata) %in% CelEsT_BM$entrezgene_id)
    WBgeneID_check <- any(row.names(DEdata) %in% CelEsT_BM$wbps_gene_id)
    genename_check <- any(row.names(DEdata) %in% CelEsT_BM$external_gene_name)
    
    shiny::validate(need(sum(c(gseq_check, entrezGeneID_check, WBgeneID_check, genename_check)) == 1,
                  "Gene IDs are inconsistent or absent; please revise. We accept WormBase sequence IDs, entrezgene IDs or WormBase IDs"))
    
    if(isFALSE(gseq_check)){
      
      if(isTRUE(entrezGeneID_check)){
        
        convertfromentrez <- CelEsT_BM[match(row.names(DEdata), CelEsT_BM$entrezgene_id), c("wormbase_gseq", "entrezgene_id")]
        uniquegseq_fromentrez <- convertfromentrez[!(duplicated(convertfromentrez$wormbase_gseq)|duplicated(convertfromentrez$wormbase_gseq, fromLast = TRUE)), "wormbase_gseq"]
        DEdata <- DEdata[match(convertfromentrez[match(uniquegseq_fromentrez, convertfromentrez$wormbase_gseq), "entrezgene_id"], row.names(DEdata)), , drop = FALSE]
        
        row.names(DEdata) <- uniquegseq_fromentrez
        
      }
      
      if(isTRUE(WBgeneID_check)){
        
        convertfromWBID <- CelEsT_BM[match(row.names(DEdata), CelEsT_BM$wbps_gene_id), c("wormbase_gseq", "wbps_gene_id")]
        uniquegseq_fromWBID <- convertfromWBID[!(duplicated(convertfromWBID$wormbase_gseq)|duplicated(convertfromWBID$wormbase_gseq, fromLast = TRUE)), "wormbase_gseq"]
        DEdata <- DEdata[match(convertfromWBID[match(uniquegseq_fromWBID, convertfromWBID$wormbase_gseq), "wbps_gene_id"], row.names(DEdata)), , drop = FALSE]
        
        row.names(DEdata) <- uniquegseq_fromWBID
        
      }
      
    } else {
      
      row.names(DEdata) <- str_remove(row.names(DEdata), "\\.[a-z]{1,2}$")
      
    }
    
    DEdata_decouple <- decoupleR::decouple(
      mat = DEdata[, "stat", drop = FALSE], 
      network = CelEsT,
      .source = "source",
      .target = "target",
      statistics = "mlm",
      args = list(mlm = list(.mor = "weight")),
      consensus_score = FALSE
    )
    
    runjs('document.getElementById("tab2status").style.display = "none";')
    
    DEdata_decouple
    
  })
  
  DEdata_forplot <- reactive({
    req(counts_computation())
    
    DEdata_decouple <- counts_computation()
    
    DEdata_decouple_wide <- pivot_wider(
      data = DEdata_decouple[, c("source", "condition", "score")], 
      names_from = "source", 
      values_from = "score"
    )
    DEdata_decouple_wide <- data.frame(t(DEdata_decouple_wide))
    
    colnames(DEdata_decouple_wide) <- DEdata_decouple_wide[1, ]
    DEdata_decouple_wide <- DEdata_decouple_wide[-1, , drop = FALSE]
    
    DEdata_decouple_wide[] <- lapply(DEdata_decouple_wide, as.numeric)
    
    DEdata_decouple_wide_p <- as.data.frame(
      pivot_wider(data = DEdata_decouple[, c("source", "condition", "p_value")], 
                  names_from = "source", values_from = "p_value")
    )
    row.names(DEdata_decouple_wide_p) <- unlist(DEdata_decouple_wide_p[, 1])
    DEdata_decouple_wide_p <- DEdata_decouple_wide_p[, -1, drop = FALSE]
    
    DEdata_decouple_wide_p <- data.frame(t(DEdata_decouple_wide_p))
    
    DEdata_decouple_wide_p[] <- lapply(DEdata_decouple_wide_p, as.numeric)
    
    DEdata_forbubbleplot <- data.frame(
      "mean_z" = unlist(DEdata_decouple_wide[row.names(DEdata_decouple_wide_p), ]),
      "p_val" = -log10(unlist(DEdata_decouple_wide_p)),
      "label" = unlist(CelEsT_BM[match(row.names(DEdata_decouple_wide_p), CelEsT_BM$wormbase_gseq), "usefullabel"])
    )
    
    DEdata_forbubbleplot[DEdata_forbubbleplot$p_val < 2.2, "label"] <- ""
    
    DEdata_forbubbleplot
    
  })
  
  TFactdata_forDT <- reactive({
    
    TFactdata_decouple <- counts_computation()
    # TFactdata_decouple <- DEdata_decouple
    
    TFactdata_decouple_wide <- pivot_wider(
      data = TFactdata_decouple[, c("source", "condition", "score")], 
      names_from = "source", 
      values_from = "score"
    )
    
    TFactdata_decouple_wide <- data.frame(t(TFactdata_decouple_wide))
    
    colnames(TFactdata_decouple_wide) <- TFactdata_decouple_wide[1, ]
    TFactdata_decouple_wide <- TFactdata_decouple_wide[-1, , drop = FALSE]
    
    TFactdata_decouple_wide[] <- lapply(TFactdata_decouple_wide, as.numeric)
    
    TFactdata_decouple_wide_p <- as.data.frame(
      pivot_wider(data = TFactdata_decouple[, c("source", "condition", "p_value")], 
                  names_from = "source", values_from = "p_value")
    )
    row.names(TFactdata_decouple_wide_p) <- unlist(TFactdata_decouple_wide_p[, 1])
    TFactdata_decouple_wide_p <- TFactdata_decouple_wide_p[, -1, drop = FALSE]
    
    TFactdata_decouple_wide_p <- data.frame(t(TFactdata_decouple_wide_p))
    
    TFactdata_decouple_wide_p[] <- lapply(TFactdata_decouple_wide_p, as.numeric)
    
    TFactdata_forDT <- data.frame(
      "Score" = unlist(TFactdata_decouple_wide[row.names(TFactdata_decouple_wide_p), ]),
      "Pval" = unlist(TFactdata_decouple_wide_p),
      "Sequence_name" = row.names(TFactdata_decouple_wide_p),
      "Gene_name" = unlist(CelEsT_BM[match(row.names(TFactdata_decouple_wide_p), CelEsT_BM$wormbase_gseq), "wormbase_locus"]),
      "Entrez_gene_id" = unlist(CelEsT_BM[match(row.names(TFactdata_decouple_wide_p), CelEsT_BM$wormbase_gseq), "entrezgene_id"]),
      "WormBase_MOR" = unlist(fullset_TFs_BM[match(row.names(TFactdata_decouple_wide_p), fullset_TFs_BM$wormbase_gseq), "Wormbase_mode_of_regulation"]),
      "UniProt_MOR" = unlist(fullset_TFs_BM[match(row.names(TFactdata_decouple_wide_p), fullset_TFs_BM$wormbase_gseq), "Uniprot_mode_of_regulation"]),
      row.names = NULL
      
    )
    
    TFactdata_forDT <- TFactdata_forDT[order(TFactdata_forDT$Pval), ]
    
    TFactdata_forDT[, c("WormBase_MOR", "UniProt_MOR")] <- str_replace_all(unlist(TFactdata_forDT[, c("WormBase_MOR", "UniProt_MOR")]), "-1", "REPRESSOR")
    TFactdata_forDT[, c("WormBase_MOR", "UniProt_MOR")] <- str_replace_all(unlist(TFactdata_forDT[, c("WormBase_MOR", "UniProt_MOR")]), "0", "BIFUNCTIONAL")
    TFactdata_forDT[, c("WormBase_MOR", "UniProt_MOR")] <- str_replace_all(unlist(TFactdata_forDT[, c("WormBase_MOR", "UniProt_MOR")]), "1", "ACTIVATOR")
    
    TFactdata_forDT
    
  })
  
  DEoutdata_forDT <- reactive({
    
    req(counts_DE())
    DEdata <- counts_DE()
    
    DEdata <- as.data.frame(DEdata)
    
    if(input$correct_speed){
      
      DEdata[, "WormBase_gene_ID"] <- row.names(DEdata)
      DEdata[, "Sequence_name"] <- CelEsT_BM[match(DEdata$WormBase_gene_ID, CelEsT_BM$wbps_gene_id), "wormbase_gseq"]
      DEdata[, "Gene_name"] <- CelEsT_BM[match(DEdata$WormBase_gene_ID, CelEsT_BM$wbps_gene_id), "wormbase_locus"]
      DEdata[, "Entrez_gene_ID"] <- CelEsT_BM[match(DEdata$WormBase_gene_ID, CelEsT_BM$wbps_gene_id), "entrezgene_id"]
      
    } else {
      
      DEdata[, "Sequence_name"] <- row.names(DEdata)
      DEdata[, "WormBase_gene_ID"] <- CelEsT_BM[match(DEdata$Sequence_name, CelEsT_BM$wormbase_gseq), "wbps_gene_id"]
      DEdata[, "Gene_name"] <- CelEsT_BM[match(DEdata$Sequence_name, CelEsT_BM$wormbase_gseq), "wormbase_locus"]
      DEdata[, "Entrez_gene_ID"] <- CelEsT_BM[match(DEdata$Sequence_name, CelEsT_BM$wormbase_gseq), "entrezgene_id"]
      
    }

    DEdata <- DEdata[, c("Sequence_name", colnames(DEdata)[1:6], "Gene_name", "WormBase_gene_ID", "Entrez_gene_ID")]

    DEdata[order(DEdata$log2FoldChange, decreasing = TRUE), ]
    
  })
  
  observe({
    if (counts_error_occurred()) {
      shinyjs::show("reset_button2")
    } else {
      hide("reset_button2")
    }
  })

  observe({
    if (input$correct_speed) {
      shinyjs::show("reset_button2")
    } else {
      hide("reset_button2")
    }
  })

  output$counts_TF_volcano <- renderPlot({
    data <- DEdata_forplot()
    
   tempcolours_for_plot <- map2color(data$p_val, pal = colorRampPalette(c("grey", "grey", "grey", "red", "red", "red"))(100))
    
    ggplot(data, aes(x = mean_z, y = p_val, size = 10^p_val, label = toupper(label))) + 
      geom_point(col = tempcolours_for_plot, alpha = 0.7) + 
      theme_classic() + 
      geom_label_repel(size = 2, max.overlaps = 40) +
      geom_vline(xintercept = 0, linetype = "dashed", col = "grey") + 
      theme(
        legend.position = "none",
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)
      ) + 
      ylab(substitute("-log"[10]~"(p-value)")) + 
      xlab("TF activity") 
  })
  
  output$counts_TF_data_table <- renderDT({
    data <- TFactdata_forDT()
    datatable(data,
              rownames = FALSE)
  }, options = list(pageLength = 10))  # You can adjust the number of rows displayed per page
  
  output$downloadcountsTFactData <- downloadHandler(
    filename = function() {
      paste("TFactivities-", Sys.Date(), ".txt", sep = "")
    },
    content = function(file) {
      write.table(TFactdata_forDT(), file, sep = "\t", row.names = FALSE, col.names = TRUE)
    }
  )
  
  output$counts_DE_data_table <- renderDT({
    req(DEoutdata_forDT())
    data <- DEoutdata_forDT()
    datatable(data,
              rownames = FALSE)
  }, options = list(pageLength = 10))  # You can adjust the number of rows displayed per page
  
  output$downloadcountsDEData <- downloadHandler(
    filename = function() {
      paste("DEdata-", Sys.Date(), ".txt", sep = "")
    },
    content = function(file) {
      write.table(DEoutdata_forDT(), file, sep = "\t", row.names = FALSE, col.names = TRUE)
    }
  )
  
  output$downloadRAPTOR <- downloadHandler(
    filename = function() {
      paste("RAPToR_age_estimates-", Sys.Date(), ".txt", sep = "")
    },
    content = function(file) {
      write.table(agecorr_process(), file, sep = "\t", row.names = FALSE, col.names = TRUE)
    }
  )
  
  observe({
    if (input$correct_speed) {
      req(agecorr_process())
      
      ageestimates_forDT <- reactive({
        
        agedata <- agecorr_process()

        agedata <- as.data.frame(agedata)
        agedata[, "sample"] <- row.names(agedata)
        agedata <- agedata[, c("sample", "age.estimate", "lb", "ub", "cor.score")]
        
        colnames(agedata) <- c("Sample", "Age estimate", "Lower bound", "Upper bound", "Corr score")
        
        agedata
        
      })
      
      
      output$age_estimates <- renderDT({
        agedata <- ageestimates_forDT()
        
        datatable(agedata, rownames = FALSE)
      }, options = list(pageLength = 10))
      
      output$age_DT_box <- renderUI({
        req(agecorr_process())
        
        shinydashboard::box(
          title = "RAPToR age estimations (hours post-hatch @ 20C)",
          status = "primary",
          solidHeader = TRUE,
          collapsible = TRUE,
          width = 12,
          downloadButton("downloadRAPTOR", "Save Data Table", class = "btn-success"),
          DTOutput("age_estimates")
        )
      })
    }
  })
  
  output$counts_volcano_ui <- renderUI({
    req(counts_computation())  # Render the UI only if the computation has been done
    shinydashboard::box(
      title = "TF Activity Volcano Plot",
      status = "primary",
      solidHeader = TRUE,
      collapsible = TRUE,
      width = 12,
      plotOutput("counts_TF_volcano", height = "500px")
    )
  })
  
  output$counts_table_ui <- renderUI({
    req(counts_computation())  # Render the UI only if the computation has been done
    shinydashboard::box(
      title = "TF activities",
      status = "primary",
      solidHeader = TRUE,
      collapsible = TRUE,
      width = 12,
      downloadButton("downloadcountsTFactData", "Save Data Table", class = "btn-success"),
      DTOutput("counts_TF_data_table")
    )
  })
  
  output$counts_DE_ui <- renderUI({
    req(counts_computation())  # Render the UI only if the computation has been done
    shinydashboard::box(
      title = "DE analysis",
      status = "primary",
      solidHeader = TRUE,
      collapsible = TRUE,
      width = 12,
      downloadButton("downloadcountsDEData", "Save Data Table", class = "btn-success"),
      DTOutput("counts_DE_data_table")
    )
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
