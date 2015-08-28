library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Welcome"),
  
  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
      #sliderInput("bins",
       #           "Number of bins:",
        #          min = 5,
         #         max = 50,
          #        value = 30)
      
      #,

    #hr(),
      #radioButtons("radio", label = h3("Radio buttons"),
      #           choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
       #          selected = 1),
    
    
    sliderInput("slider2", label = h4("Range for sperm"), min = 0, 
                max = 10, value = c(8, 10)),
    
      sliderInput("slider3", label = h4("Range for spore samples"), min = 0, 
                max = 10, value = c(0, 10)),
    
    sliderInput("slider1", label = h4("Range for all but sperm and spore samples"), min = 0, 
                max = 10, value = c(0, 5)),
    
   checkboxInput("outliers", label="Exclude variance outliers",value= FALSE)),
   
   #hr(),
   #fluidRow(column(3, verbatimTextOutput("value"))),
  
    # Show a plot of the generated distribution
    mainPanel(
      #plotOutput("distPlot"),
      #plotOutput("bxPlt"),
      plotOutput("matplot"),
      dataTableOutput('mytable')
    )
  )
))
