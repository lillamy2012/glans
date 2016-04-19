library(shiny)
library(DESeq2)

load("TSSdds.Rdata")

tpm=read.csv("abundance_per_tx.csv",row.names=1)
res = results(dds)
resOrdered <- res[order(res$padj<0.1,abs(res$"log2FoldChange"),decreasing=T),]
counts = counts(dds,normalized=T)
colnames(counts) = gsub("/","", sapply(strsplit(colnames(counts),"\\."),"[[",3))
scn = asinh(rowMeans(counts[,1:2]))
vcn =  asinh(rowMeans(counts[,3:4]))
sig = which(res$padj<0.1)




shinyServer(function(input, output) {


 # output$scatterATAC <- renderPlot({
#    plot(scn,vcn,pch=".")
 #   points(scn[sig],vcn[sig],col="red",pch=".")
#    abline(0,1,col="red")
 #   })
  
  output$barplotATAC = renderPlot({
    s = input$table_rows_selected
    print(s)
    if (length(s)){
      barplot(counts[s,],las=2,ylab="norm.counts")
    } 
  })
  
  output$scatterATACsel = renderPlot({
    s = input$table_rows_selected
    print(s)
    plot(scn,vcn,pch=".",col="grey",ylab="asinh(vsn)",xlab="asinh(scn)")
    points(scn[sig],vcn[sig],col="red",pch=".")
    if (length(s)){
      points(scn[s],vcn[s],col="blue")
    }
  })
  
  output$histRNAsel = renderPlot({
    s = input$table_rows_selected
    print(s)
    hist(asinh(tpm$amean),n=20,main="",ylab="",xlab="asinh(tpm)")
    if (length(s)){
      for (i in s){
        abline(v=asinh(tpm[s,"amean"]),col="red")
      }
    }
  })
  
  output$barplotRNA = renderPlot({
    s = input$table_rows_selected
    print(s)
    if (length(s)){
      barplot(as.matrix(tpm[s,1:3]),las=2)
      print(as.matrix(tpm[s,1:3]))
    }
  })
  
  
  
output$table <- DT::renderDataTable(
  as.data.frame(resOrdered)
)


#output$info <- renderText({
#  paste0("x=", input$plot_click$x, "\ny=", input$plot_click$y)
#})
}

)