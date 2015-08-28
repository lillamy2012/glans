library(shiny)
library(DT)

server <- function(input, output) {
  output$iris_type <- DT::renderDataTable({
    datatable(data.frame(Species=paste0("<a href='#filtered_data'>", unique(iris$Species), "</a>")),
              escape = FALSE,
              callback = JS(
                'table.on("click.dt", "tr", function() {
                tabs = $(".tabbable .nav.nav-tabs li a");
                $(tabs[1]).click();})'))
})
  
  output$filtered_data <- DT::renderDataTable({
    selected <- input$iris_type_rows_selected
    if(is.null(selected)){
      datatable(iris)
    } else {
      datatable(iris[iris$Species %in% unique(iris$Species)[selected], ])
    }
  })
}

ui <- shinyUI(fluidPage(
  mainPanel(
    tabsetPanel(
      tabPanel("Iris Type", DT::dataTableOutput("iris_type")),
      tabPanel("Filtered Data", DT::dataTableOutput("filtered_data"))
    )
  )
))

shinyApp(ui = ui, server = server)