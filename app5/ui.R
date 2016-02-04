
#library(DESeq2)
library(shiny)

shinyUI(fluidPage(
  titlePanel("Mp groups"),
  sidebarLayout(
    sidebarPanel(
      checkboxGroupInput("outlier", "Remove",
                         c( "young_male_reprod_r1" = "young_male_reprod_r1",
                            "young_male_reprod_r2" = "young_male_reprod_r2",
                            "Mp_WT_sperm" = "Mp_WT_sperm",
                            "Mp_sperm1" = "Mp_sperm1",
                            "Mp_sperm2" ="Mp_sperm2",
                            "Mp_sperm3" = "Mp_sperm3", 
                            "10DAF_r1" = "10DAF_r1",
                            "10DAF_r2" ="10DAF_r2",
                            "15DAFembry1" = "15DAFembry1",
                            "15DAFembry2" = "15DAFembry2",
                            "20DAF1mmembry1" = "20DAF1mmembry1",
                            "20DAF1mmembry2" = "20DAF1mmembry2",
                            "20DAFSoftembry1"= "20DAFSoftembry1",
                            "20DAFSoftembry2" = "20DAFSoftembry2",
                            "spore1" =  "spore1",
                            "spore2" = "spore2",
                            "female_3mm_r1" = "female_3mm_r1",
                            "female_3mm_r2" = "female_3mm_r2",
                            "Mp_female" = "Mp_female",     
                            "Male_thalli_r1" = "Male_thalli_r1", 
                            "Male_thalli_r2" =   "Male_thalli_r2",
                            "Mp_thalli" = "Mp_thalli"
                              ),
                         
                         selected = NULL, inline = FALSE, width = NULL),
      actionButton("goButton", "Go!"),
      radioButtons("Group","Group",
                   c("Sperm/young rep" = 1,
                   "15DAFembryo/20DAF1mm" = 2,
                   "20DAF1mm/20DAFsoft" = 3,
                   "20DAFsoft/spore" =4,
                   "10DAF/15DAFembryo"=5)),
      
      radioButtons("Set","Set",
                   c("1 or 2" = "union",
                     "1 and 2" = "intersect",
                     "1 and not 2" = "setdiff1",
                     "2 and not 1" = "setdiff2")),
    
     numericInput("padj", label = h4("Max Anova p-value (adjusted)"), value = 0.1,max=1,min=0,step=0.001),
     numericInput("minFC", label = h4("Min fold change"), value = 2),
     numericInput("minavFC", label = h4("Min average fold change"), value = 2),
     checkboxInput("transf",label="asinh transform data",value=TRUE)
    # textInput("collection_txt",label="Foo")
     ),

    mainPanel(plotOutput("matplot",click="plot_click"),verbatimTextOutput("click_info"),
              hr(),
              #textOutput("chos"),
              plotOutput("selc"),
              downloadButton('downloadData', 'Download list as it is'),
              dataTableOutput("tot_data")
             ))))
