library(shiny)
library(dplyr)
library(DT)
source("functions.R")
library(ggplot2)
library(stringr)
###################################
shinyServer(function(input, output) {

  readIn1 <- reactive({
  inFile <- input$file1
    if (is.null(inFile))
    return(NULL)
  dd = read.csv(inFile$datapath, header=TRUE, sep=input$sep, 
                quote=input$quote,skip=2,dec=",")
  #dd = convert(dd)
  return(dd)
  })
  
  readIn2 <- reactive({
    inFile <- input$file2
    if (is.null(inFile))
      return(NULL)
    dd = read.csv(inFile$datapath, header=TRUE, sep=input$sep, 
                  quote=input$quote,skip=2,dec=",",comment="")
    extra = dd[1,]
    dd = dd[-1,]
    dd = dd[!is.na(dd[,1]),]
    return(dd)
  })
  
  
  output$summary =  DT::renderDataTable({
    dd= readIn2()
    indata=readIn1()
    indata=filter_amanda(indata,input$filter)
    dd= fasta_format(dd)
    res = matrix(NA,nrow=nrow(dd),ncol=4)
    withProgress(message = 'Calculating', value = 0, {
      for (i in 1:nrow(dd)){
      fasta=dd$Sequence[i]
      index=grep(dd[i,"Accession"],indata$Accession,fixed=T)
      indind=indata[index,]
      res[i,]=unlist(summaryFunction(ProteinPlotMat(fasta,indind)[[1]],indind,dd[i,"Accession"]))
      # Increment the progress bar, and update the detail text.
      incProgress(1/nrow(dd), detail = paste("Protein", dd[i,"Accession"]))
      
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
      }
})
    as.data.frame(res)
    rownames(res)=dd$Accession
    colnames(res) = names(summaryFunction(ProteinPlotMat(fasta,indind)[[1]],indind,dd[i,"Accession"]))
    res
  },options = list(lengthMenu = c(5, 30, 50), pageLength = 15),escape=FALSE,callback = JS(
    'table.on("click.dt", "tr", function() {
    tabs = $(".tabbable .nav.nav-tabs li a");
    $(tabs[1]).click();})'))
  
  output$ProteinList = DT::renderDataTable({
    dd= readIn1()
    if (is.null(dd))
      return(NULL)
    gr = sapply(strsplit(as.character(dd$Accession),";"),length)
    dd$unq = ifelse(gr>1,0,1)
    ll = data.frame(Protein=unique(dd$Accession[which(dd$unq==1)]))
    ll
  })
  
  output$FastaList = DT::renderDataTable({
    dd= readIn2()
    if (is.null(dd))
      return(NULL)
    dd = fasta_format(dd)
    dd$Accession = paste0("<a href='#filtered_data'>", dd$Accession, "</a>")
    dd
  },options = list(lengthMenu = c(5, 30, 50), pageLength = 15),escape=FALSE,callback = JS(
    'table.on("click.dt", "tr", function() {
    tabs = $(".tabbable .nav.nav-tabs li a");
    $(tabs[1]).click();})'))
  
  output$my_prot = renderPlot({
    dd= readIn2()
    indata=readIn1()
    indata=filter_amanda(indata,input$filter)
    dd= fasta_format(dd)
    selected=input$FastaList_rows_selected
    selected = selected[length(selected)]
    index=grep(dd[(selected),"Accession"],indata$Accession,fixed=T)
    indind=indata[index,]
    if(nrow(indind)>0){
      mm = ProteinPlotMat(dd[as.numeric(selected),"Sequence"],indind)
      
      ProteinPlot(mm[[1]],mm[[2]])
    }
    else
      plot(1:100,pch=NULL)
      })
  
  output$info = renderDataTable({
    dd= readIn2()
    indata=readIn1()
    dd= fasta_format(dd)
    selected=input$FastaList_rows_selected
    selected = selected[length(selected)]
    index=grep(dd[(selected),"Accession"],indata$Accession,fixed=T)
    indind=indata[index,]
    mm = ProteinPlotMat(dd[as.numeric(selected),"Sequence"],indind)
    info = summaryFunction(mm[[1]],indind,dd[(selected),"Accession"])
    return(info)
    },options = list(searching = FALSE,paging = FALSE))
    
  output$mod = renderDataTable({
    dd= readIn2()
    indata=readIn1()
    dd= fasta_format(dd)
    selected=input$FastaList_rows_selected
    selected = selected[length(selected)]
    index=grep(dd[(selected),"Accession"],indata$Accession,fixed=T)
    indind=indata[index,]
    mm = ProteinPlotMat(dd[as.numeric(selected),"Sequence"],indind)
    res=(modSummary(mm[[1]]))
    colnames(res) = c("Position","Modification","# Modifications","Procentage of all peptides")
    return(res)
  })
  
})