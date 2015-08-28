library(DT)
shinyUI({
  basicPage(

    
    tags$div(id="blah", class="shiny-input-radiogroup", 
             dataTableOutput( "results" )
    ), 
    tags$h1(textOutput( "foo" ))
  )
})
