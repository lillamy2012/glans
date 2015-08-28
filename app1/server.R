library(shiny)
load("../data/app1_data.rdata")
outl=(unique(which(log10(sMat[,-c(1, 6, 9, 11)])>1,arr.ind=T))[,1])
shinyServer(function(input, output) {
  
  dataInput <- reactive({
    if (is.null(input$slider1[2]) || is.na(input$slider1[1])){
      return()
    } else {
      mostH = input$slider1[2]
      mostL = input$slider1[1]
      spermH = isolate(input$slider2[2])
      spermL = isolate(input$slider2[1])
      sporL= isolate(input$slider3[1])
      sporH = isolate(input$slider3[2])
      index=which(qdata[,1]<=mostH & qdata[,1]>mostL & 
                    qdata[,2]<=mostH & qdata[,2]>=mostL &
                    qdata[,3]<=mostH & qdata[,3]>=mostL &
                    qdata[,4]<=mostH & qdata[,4]>=mostL &
                    qdata[,5]<=mostH & qdata[,5]>=mostL &
                    qdata[,6]<=mostH & qdata[,6]>=mostL &
                    qdata[,7]<=mostH & qdata[,7]>=mostL &
                    qdata[,8]<=mostH & qdata[,8]>=mostL &
                    qdata[,9]<=mostH & qdata[,9]>=mostL &
                    qdata[,10]<=mostH & qdata[,10]>=mostL &
                    qdata[,11]<=mostH & qdata[,11]>=mostL &
                    qdata[,12]<=mostH & qdata[,12]>=mostL &
                    qdata[,13]<=mostH & qdata[,13]>=mostL &
                    qdata[,14]<=mostH & qdata[,14]>=mostL &
                    qdata[,15]<=spermH & qdata[,15]>=spermL &
                    qdata[,16]<=sporH & qdata[,16]>=sporL &
                    qdata[,17]<=sporH & qdata[,17]>=sporL & 
                    qdata[,18]<=mostH & qdata[,18]>=mostL & 
                    qdata[,19]<=mostH & qdata[,19]>=mostL)
    return(index)
  }})

                 
  output$matplot <-renderPlot({
    index=dataInput()
    x    <- as.matrix(asinh(meansMat[,-1]))
    #print(index)
    matplot(t(x[index,]),type="b",pch=19,lty=1)
    
  })
  output$value <- renderPrint({ input$checkbox })
  output$mytable = renderDataTable({
    index=dataInput()
    if(input$outliers){
        print("ok")
      (meansMat[setdiff(index,outl),])
    } else 
    #index=dataInput()
    #print(summary(as.matrix((meansMat[,-1]))))
    #print(setdiff(index,input$outliers))
    (meansMat[index,])
  
  })
})