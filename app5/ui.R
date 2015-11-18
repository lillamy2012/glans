library(shiny)

shinyUI(fluidPage(
  titlePanel("Mp groups"),
  sidebarLayout(
    sidebarPanel(
      radioButtons("Group","Group",
                   c("Sperm/young rep" = 1,
                   "15DAFembryo/20DAF1mm" = 2,
                   "20DAF1mm/20DAFsoft" = 3,
                   "20DAFsoft/spore" =4)),
      
      radioButtons("Set","Set",
                   c("1 or 2" = "union",
                     "1 and 2" = "intersect",
                     "1 and not 2" = "setdiff1",
                     "2 and not 1" = "setdiff2")),
    
     numericInput("padj", label = h4("Max Anova p-value (adjusted)"), value = 0.1,max=1,min=0,step=0.001),
     numericInput("minFC", label = h4("Min fold change"), value = 2),
     numericInput("minavFC", label = h4("Min average fold change"), value = 2),
     checkboxInput("transf",label="asinh transform data",value=TRUE)
     ),

    mainPanel(plotOutput("matplot",click="plot_click"),verbatimTextOutput("click_info"),
              hr(),
              plotOutput("selc"),
              downloadButton('downloadData', 'Download list as it is'),
              dataTableOutput("tot_data")
             ))))