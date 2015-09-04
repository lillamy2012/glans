library(shiny)

shinyUI(fluidPage(
  titlePanel("Protein analysis app (PaApp)"),
  sidebarLayout(
    sidebarPanel(
      fileInput('file1', 'Choose CSV File',
                accept=c('text/csv', 
                         'text/comma-separated-values,text/plain', 
                         '.csv')),
      #checkboxInput('header', 'Header', TRUE),
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
                  max = 1000, value = 200),
      tags$hr(),
        fileInput('file2', 'Choose CSV File for FASTA',
                  accept=c('text/csv', 
                           'text/comma-separated-values,text/plain', 
                           '.csv'))
        #checkboxInput('header', 'Header', TRUE),
       # radioButtons('sep', 'Separator',
        #             c(Comma=',',
         #              Semicolon=';',
        #               Tab='\t'),
         #            ';'),
        #radioButtons('quote', 'Quote',
         #            c(None='',
          #             'Double Quote'='"',
           #            'Single Quote'="'"),
            #         '"')
       
    ),
  mainPanel(
      tabsetPanel(
        tabPanel("Fasta List", DT::dataTableOutput("FastaList")),
        tabPanel("Plot Data", plotOutput("my_prot"),dataTableOutput("info"),dataTableOutput("mod")),
        tabPanel("Summary", dataTableOutput("summary"))
      )
  )
    )))

