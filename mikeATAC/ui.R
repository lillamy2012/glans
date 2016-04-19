library(shiny)

shinyUI(fluidPage(
  titlePanel("ATAC-seq"),
    #plotOutput(outputId = "scatterATAC", click="plot_click"),verbatimTextOutput("info"),
    fluidRow(
    column(width = 3, class = "well",
           h4("ATAC-seq counts"),
    plotOutput(outputId = "scatterATACsel")),
    column(width = 3, class = "well",
           h4("ATAC-seq per replicate (counts)"),
          plotOutput(outputId = "barplotATAC")),
    column(width = 3, class = "well",
           h4("RNA-seq TPM"),
           plotOutput(outputId = "histRNAsel")),
    column(width = 3, class = "well",
           h4("RNA-seq per replicate (TPM)"),
           plotOutput(outputId = "barplotRNA"))
    
    ),
    
    DT::dataTableOutput("table")
))
    
  