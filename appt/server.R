library(DT)
server.fun <- function(input, output) {
  output$foo <- renderText({
    input$blah
  })
  
  output$results <- renderDataTable({
   
      results <- data.frame(ID=1:10, contents=letters[1:10])
      results$checkbox <- sprintf( '<input type="radio" name="blah" value="%d"/>', 1:10 )
      datatable(results, escape=FALSE)
  
    
  })
}
