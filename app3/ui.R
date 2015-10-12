library(shiny)
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
      
      conditionalPanel(
        condition = "output.done",
      checkboxInput("grouping", label = "Combine all samples", value = TRUE)),
      
    
      conditionalPanel(
        condition = "output.group",
        uiOutput("group1")),
            
      conditionalPanel(
        condition = "output.group",
        uiOutput("group2"))),
    
  mainPanel(
      tabsetPanel(
        tabPanel("Fasta List", DT::dataTableOutput("FastaList")),
        tabPanel("", plotOutput("my_prot"),dataTableOutput("info"),hr(),hr(),dataTableOutput("mod"),hr(),downloadButton('downloadData', 'Download full list')),
        tabPanel("Summary", dataTableOutput("summary")),
        tabPanel("", plotOutput("my_prot1"),dataTableOutput("info1"),hr(),hr(),dataTableOutput("mod1"),hr(),downloadButton('downloadData1', 'Download full list'))
      
       )
  )
    )))

