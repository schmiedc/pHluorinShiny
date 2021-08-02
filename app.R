packages <- c("gridExtra", "shiny", "shinyFiles", "reshape2", "plyr", "tidyverse", "broom", "ggplot2", "rstudioapi")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

library(gridExtra)
library(shiny)
library(shinyFiles)

# specify the location where to find scripts
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("dataProcessing.R")
source("saveData.R")
source("plotData.R")
source("fitting.R")
# ============================================================================
#
#  DESCRIPTION: Data analysis for Fíji pHlorin workflow
#              
#       AUTHOR: Christopher Schmied, 
#      CONTACT: schmied@dzne.de
#     INSITUTE: Leibniz-Forschungsinstitut f r Molekulare Pharmakologie (FMP)
#               Cellular Imaging - Core facility
#               Campus Berlin-Buch
#               Robert-Roessle-Str. 10
#               13125 Berlin, Germany
#
#         BUGS:
#        NOTES: 
# DEPENDENCIES: ggplot2: install.packages("ggplot2")
#               reshape2: install.packages("reshape2")
#               plyr: install.packages("plyr")
#               gridExtra: install.packages("gridExtra")
#               xlsx: install.packages("xlsx")
#
#      VERSION: 1.0.0
#      CREATED: 2018-05-24
#     REVISION: 2018-08-07
#
# ============================================================================

ui <- fluidPage( 
  
  titlePanel("pHluorin Data Processing"),
  
  sidebarLayout(position = "right",
  
    sidebarPanel(
      
      fluidRow(
        
        h3("Settings:"),
        
        shinyDirButton("inputdir", "Input directory", "Please select a folder"),
        verbatimTextOutput("directorypath"),
        
        numericInput(inputId = "timeRes", label = "Set frame rate:", 2, min = 0, max = 10, step =01,
                     width = NULL),
        
        numericInput(inputId = "frameStim", label = "Stimulation Frame:", 5, min = 0, max = 100, step =01,
                     width = NULL),
        
        textInput(inputId = "resultname", label = "Result Name:", value = "", width = NULL,
                  placeholder = "Test")

      ),
      
      fluidRow(
        
        h3("Processing:"),
        
        actionButton("processData", "Load Data"),
        
        h3("Save Results:"),
        
        actionButton("saveData", "Save Results")
        #checkboxInput(inputId = "asCsv", "as .csv", value = FALSE, width = NULL),
        #checkboxInput(inputId = "asXlsx", "as .xlsx", value = FALSE, width = NULL)
        
      ),
      
      fluidRow(
        
        tags$hr(),
        tags$p("Author: Christopher Schmied"),
        tags$p("Contact: schmied@fmp-berlin.de"),
        tags$p("Cellular Imaging - Core facility"),
        tags$hr()
        
      ),
      
      fluidRow(
      tags$img(width = 450, src='FMP_Wort-Bildmarke_nebeneinander_gruen_RGB.png')
      )
      
      ),
    
    mainPanel(
      
      tabsetPanel(type = "tabs",
      
        tabPanel("Overview", plotOutput("summary2", height=1200)),
        tabPanel("Detail", selectInput(inputId = "dataName", label = "Input Name:", 
                                       choices = "", selected = ""),
                 actionButton("plotDetail", "plot Detail"),
                 plotOutput("detail", height=1200)
                 )
        
      )
      
    )
  
  )
  
)

server <- function(input, output, session) {
  
  volumes <- c(Home = fs::path_home(), 
               "R Installation" = R.home(), 
               getVolumes()())
  
  # dir
  shinyDirChoose(
    input, 
    'inputdir', 
    roots = volumes,
    session = session, 
    restrictions = system.file(package = "base"))
  
  global <- reactiveValues(datapath = getwd())

  dir <- reactive(input$inputdir)
  global <- reactiveValues(datapath = getwd())
  
  observe({
    cat("\ninput$directory value:\n\n")
    print(input$inputdir)
  })
  
  output$directorypath <- renderPrint({
    
    if (is.integer(input$inputdir)) {
      cat("No directory has been selected (shinyDirChoose)")
    } else {
      global$datapath <- parseDirPath(volumes, input$inputdir)
      parseDirPath(volumes, input$inputdir)

    }
  })
  
  observeEvent(input$plotDetail, {
    
    tryCatch({
      
      # further settings
      labelSignal = "Spot"
      labelBackground = "background"
      
      inputDirectory <- global$datapath
      
      # collects raw data
      table.signal <- collectList(inputDirectory, labelSignal, input$timeRes)
      table.background <- collectList(inputDirectory, labelBackground, input$timeRes)
      
      detail <- list()
      detail[["one"]] <- plotRawMean_single(table.signal, input$dataName)
      detail[["two"]] <- plotRawMean_single(table.background, input$dataName)
      
      # generate output plots
      output$detail <- renderPlot({
        
        plots <- grid.arrange(grobs = detail, ncol = 1, top = "Raw traces")
        
        print(plots)
        
      })
      
    }, error=function(e) {
      
      message(e)
      showNotification(paste0("WARNING:   ", e), type = 'error')
      
    }, warning=function(w) {
      
      message(w)
      showNotification(paste0("WARNING:   ", w), type = 'warning')
      
    })
    
})
  
  # script
  # gets lists for signal and background
  observeEvent(input$processData, {
    
    tryCatch({
    
      # further settings
      labelSignal = "Spot"
      labelBackground = "background"
      
      inputDirectory <- global$datapath
      
      # collects raw data
      table.signal <- collectList(inputDirectory, labelSignal, input$timeRes)
      table.background <- collectList(inputDirectory, labelBackground, input$timeRes)
      
      updateSelectInput(session, 
                        "dataName", 
                        label = NULL, 
                        choices = unique(table.signal$name),
                        selected = NULL
                        )
      
      # calculates average mean intensity per frame
      avg.signal <- calcMean(table.signal)
      avg.background <- calcMean(table.background)

      finalTable <- processData(indir, input$frameStim, avg.signal, avg.background)
      
      summary2 <- plotSummary2(avg.signal, avg.background, finalTable)
      
      output$summary2 <- renderPlot({
        
        plots <- grid.arrange(grobs = summary2, ncol = 1, top = "Raw Mean per Movie")
        
        print(plots)
        
      })
      
      
    }, error=function(e) {
      
      message(e)
      showNotification(paste0("WARNING:   ", e), type = 'error')
      
    }, warning=function(w) {
      
      message(w)
      showNotification(paste0("WARNING:   ", w), type = 'warning')
      
    })
      
  })
    
  observeEvent(input$saveData, {
      
      tryCatch({
      
        inputDirectory <- global$datapath
        outputDirectory <- global$datapath
        resultname <- input$resultname
        
        # further settings
        labelSignal = "Spot"
        labelBackground = "background"
        
        inputDirectory <- global$datapath
        
        # collects raw data
        table.signal <- collectList(inputDirectory, labelSignal, input$timeRes)
        table.background <- collectList(inputDirectory, labelBackground, input$timeRes)
        
        # calculates average mean intensity per frame
        avg.signal <- calcMean(table.signal)
        avg.background <- calcMean(table.background)
        
        # generate final table
        finalTable <- processData(indir, input$frameStim, avg.signal, avg.background)
        
        tau <- calcTau(finalTable)
        
        # save files
        writeToCsv(outputDirectory, resultname, table.signal, table.background, finalTable, tau)
        
        # plot data
        # ======================================================================
        # create result plots
        raw_signal <- plotRawMean(table.signal)
        raw_signal_grids <-  marrangeGrob(raw_signal, ncol = 3, nrow = 4, top = "Raw grey values of active boutons")
        ggsave(plot = raw_signal_grids,
               file=file.path(outputDirectory, paste0(resultname, "_rawPlotsSignal.pdf") ), 
               width = 297, 
               height = 210, 
               units = "mm") 
        
        # plots Raw grey values and area of background
        raw_signal <- plotRawMean(table.background)
        raw_signal_grids <-  marrangeGrob(raw_signal, ncol = 3, nrow = 4, top = "Raw grey values of background")
        ggsave(plot = raw_signal_grids,
               file=file.path(outputDirectory, paste0(resultname, "_rawPlotsBackground.pdf") ), 
               width = 297, 
               height = 210, 
               units = "mm") 
        
        # RawAreaBoutons
        rawArea <- plotRawArea(table.signal, "active boutons")
        ggsave(plot = rawArea,
               file=file.path(outputDirectory, paste0(resultname, "_rawAreaBoutons.pdf") ),
               width = 297, 
               height = 210, 
               units = "mm") 
        
        # RawAreaBackground
        rawAreaBackground <- plotRawArea(table.background, "background")
        ggsave(plot = rawAreaBackground,
               file=file.path(outputDirectory, paste0(resultname, "_rawAreaBackground.pdf") ),
               width = 297, 
               height = 210, 
               units = "mm")  
        
        # Plots mean values per movie
        plot.list <- plotMeans(avg.signal, avg.background, finalTable)
        
        test_plots <- marrangeGrob(plot.list, ncol = 1, nrow = 1, top = "Processing results")
        
        ggsave(plot = test_plots,
               file=file.path(outputDirectory, paste0(resultname, "_meanResults.pdf") ),
               width = 297, 
               height = 210, 
               units = "mm") 

        showNotification("Data saved")
        
      }, error=function(e) {
        
        message(e)
        showNotification(paste0("WARNING:   ", e), type = 'error')
        
      }, warning=function(w) {
        
        message(w)
        showNotification(paste0("WARNING:   ", w), type = 'warning')
        
      })
    
  })
  
}

# Run the application
shinyApp(ui = ui, server = server)