library(shiny)
library(DT)
shinyUI(fluidPage(
  titlePanel("Protein analysis app (PaApp)"),
  sidebarLayout(
    sidebarPanel(
      fileInput('file1', 'Choose CSV File with PSM details',
                accept=c('text/csv', 
                         'text/comma-separated-values,text/plain', 
                         '.csv')),
      tags$hr(),
      fileInput('file2', 'Choose CSV File with Coverage details',
                accept=c('text/csv', 
                         'text/comma-separated-values,text/plain', 
                         '.csv')),
  
      numericInput("filter", label = h4("Min. Amanda score"), value = 100),
      hr(),
      checkboxInput("checkbox", label = "Use ONLY peptides UNIQUE mapped to protein", value = FALSE),
      checkboxInput("norm", label = "Normalize based on total number peptides per sample", value = TRUE),
      
      conditionalPanel(
        condition = "output.done",
      checkboxInput("grouping", label = "Combine all samples", value = TRUE)),
      
    
      conditionalPanel(
        condition = "output.group",
        uiOutput("group1")),
            
      conditionalPanel(
        condition = "output.group",
        uiOutput("group2")),
      
      conditionalPanel(
        condition = "output.group",
        actionButton("goButton", "Go!"))),
    
  mainPanel(
      tabsetPanel(
        tabPanel("Fasta List", DT::dataTableOutput("FastaList")),
        
        tabPanel("", #conditionalPanel(
                #condition = "output.track",
                 plotOutput("my_prot",height="600px")),
                #dataTableOutput("info"),
                 #hr(),hr(),dataTableOutput("mod"),hr(),
                 #downloadButton('downloadData', 'Download full list')),
        
        tabPanel("Summary", dataTableOutput("summary"),
                 downloadButton('downloadSummary', 'Download full list')),
        
        tabPanel("", plotOutput("my_prot1",height="600px"),
                 checkboxInput('returnpdf', 'output pdf', FALSE),
                 conditionalPanel(
                   condition = "input.returnpdf == true",
                   fluidRow(
                     column(3,
                            downloadButton('downloadData0', 'Download figure')
                     ),
                   column(4,
                          wellPanel(
                          sliderInput('height', 'Select pdf height', 3, 10, 7))),
                   column(4,
                          wellPanel(
                          sliderInput('size', 'Select pdf width', 6, 25, 15))
                 ))),
                 dataTableOutput("info1"),
                 hr(),hr(),dataTableOutput("mod1"),hr(),
                 downloadButton('downloadData1', 'Download full list'),
                 dataTableOutput("mod2"),hr(),
                 conditionalPanel(
                  condition = "output.group",
                  downloadButton('downloadData2', 'Download full list'))),
        
        tabPanel("ScatterPlot",plotOutput("plotGroups",height="600px",
                width = "600px",click="plot_click"),verbatimTextOutput("click_info"))
       )
  )
    )))

