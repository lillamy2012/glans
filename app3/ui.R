library(shiny)

shinyUI(fluidPage(
  
  # Application title
  titlePanel("Protein analysis app (PaApp)"),
  
  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
      fileInput('file1', 'Choose CSV File',
                accept=c('text/csv', 
                         'text/comma-separated-values,text/plain', 
                         '.csv')),
      
      
      checkboxInput('header', 'Header', TRUE),
      radioButtons('sep', 'Separator',
                   c(Comma=',',
                     Semicolon=';',
                     Tab='\t'),
                   ';'),
      radioButtons('quote', 'Quote',
                   c(None='',
                     'Double Quote'='"',
                     'Single Quote'="'"),
                   '"'),
      tags$hr(),
      sliderInput("filter", label = h4("Amanda score filter"), min = 100, 
                  max = 1000, value = 200)),
     
  mainPanel(
      tabsetPanel(
        tabPanel('Protein List',DT::dataTableOutput("test")),
        tabPanel("Scatter plots",plotOutput("plot",click="plot_click"),
                 plotOutput("plot1",click="plot_click"),
                 verbatimTextOutput("click_info")),
        tabPanel("Histograms",plotOutput("hist")),
        tabPanel("Group", DT::dataTableOutput("mytable")),
        tabPanel("Filtered Data", DT::dataTableOutput("filtered_data"))
      )
  )
)))
