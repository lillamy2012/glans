library(shiny)
library(dplyr)
library(DT)
source("functions.R")
library(stringr)
###################################
out=NULL
shinyServer(function(input, output) {

################################################################
### read in files functions 1 = psm details, 2 = coverage details  
################################################################
  
  readIn1 <- reactive({
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    dd = read.csv(inFile$datapath, header=TRUE, sep=";", 
                quote='"',skip=2,dec=",",stringsAsFactors=FALSE)
    dd = convert(dd)
    return(dd)
  })
  
  readIn2 <- reactive({
    inFile <- input$file2
    if (is.null(inFile))
      return(NULL)
    dd = read.csv(inFile$datapath, header=TRUE, sep=";", 
                  quote="'",skip=2,dec=",",comment="",stringsAsFactors=FALSE)
    #extra = dd[1,]
    dd = dd[-1,]
    dd = dd[!is.na(dd[,1]),]
    return(dd)
  })

##################################################################    
## summarize peptides, modifications etc per protein (tab 'Summary')
##################################################################  
  
  output$summary =  DT::renderDataTable({
    dd= readIn2()
    indata=readIn1()
    if(is.null(indata)|is.null(dd))
      return(NULL)
    indata=filter_amanda(indata,input$filter)
    dd= fasta_format(dd)
    res = matrix(NA,nrow=nrow(dd),ncol=4)
    withProgress(message = 'Calculating', value = 0, {
      for (i in 1:nrow(dd)){
      fasta=dd$Sequence[i]
      if(input$checkbox==TRUE){
        selected=sub("\\|","\\\\|",dd[i,"Accession"])
        grepWord=paste("^","$",sep=selected)
        index = grep(grepWord,indata$Accession,fixed=F,value=F)
      } else {
        index=grep(dd[i,"Accession"],indata$Accession,fixed=T)
      }
      indind=indata[index,]
      if(length(index)>0)
        res[i,]=unlist(summaryFunction(ProteinPlotMat(fasta,indind)[[1]],indind,dd[i,"Accession"]))
      # Increment the progress bar, and update the detail text.
      incProgress(1/nrow(dd), detail = paste("Protein", dd[i,"Accession"]))
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
      }
})
    as.data.frame(res)
    rownames(res)=dd$Accession
    #colnames(res) = names(summaryFunction(ProteinPlotMat(fasta,indind)[[1]],indind,dd[i,"Accession"]))
    colnames(res)=c("coverage","nr.peptides","modification_types","modifications_detected")
    res
  },options = list(lengthMenu = c(5, 30, 50), pageLength = 15),escape=FALSE,callback = JS(
    'table.on("click.dt", "tr", function() {
    tabs = $(".tabbable .nav.nav-tabs li a");
    $(tabs[3]).click();})'))

  ##################################################################    
  ## summarize proteins in study (tab 'Fasta List')
  ##################################################################  
  
  output$FastaList = DT::renderDataTable({
    dd= readIn2()
    if (is.null(dd))
      return(NULL)
    dd = fasta_format(dd)
    #dd$Accession = paste0("<a href='#filtered_data'>", dd$Accession, "</a>")
    dd
  },options = list(lengthMenu = c(5, 30, 50), pageLength = 15),escape=FALSE,callback = JS(
    'table.on("click.dt", "tr", function() {
    tabs = $(".tabbable .nav.nav-tabs li a");
    $(tabs[1]).click();})'))
  
  ##################################################################    
  ## plot protein from Fasta List (tab 1)
  ##################################################################  
  
  output$my_prot = renderPlot({
    dd= readIn2()
    indata=readIn1()
    indata=filter_amanda(indata,input$filter)
    dd= fasta_format(dd)
    selected=input$FastaList_rows_selected
    selected = selected[length(selected)]
    if(input$checkbox==TRUE){
      index = use_unique(dd,selected,indata)
    } else {
    index=grep(dd[(selected),"Accession"],indata$Accession,fixed=T)
    }
    indind=indata[index,]
    if(nrow(indind)>0){
      mm = ProteinPlotMat(dd[as.numeric(selected),"Sequence"],indind)
      
      ProteinPlot(mm[[1]],mm[[2]])
    }
    else
      plot(1:100,pch=NULL)
      })
  ##################################################################    
  ## plot protein from Summary (tab 3)
  ##################################################################  
  
  output$my_prot1 =renderPlot({
    dd= readIn2()
    indata=readIn1()
    indata=filter_amanda(indata,input$filter)
    dd= fasta_format(dd)
    selected=input$summary_rows_selected
    selected = selected[length(selected)]
    selected=sub("\\|","\\\\|",selected)
    sel=(grep(selected,dd[,"Accession"]))
    if(input$checkbox==TRUE){
      index = use_unique(dd,sel,indata)
    } else {
      index=grep(dd[(sel),"Accession"],indata$Accession,fixed=T)
    }
    indind=indata[index,]
    if(nrow(indind)>0){
      mm = ProteinPlotMat(dd[as.numeric(sel),"Sequence"],indind)
      ProteinPlot(mm[[1]],mm[[2]])
    }
    else
      plot(1:100,pch=NULL)
  })
    
  ##################################################################    
  ## info summary from Fasta List (tab 1)
  ##################################################################  
  
  output$info = renderDataTable({
    dd= readIn2()
    indata=readIn1()
    indata=filter_amanda(indata,input$filter)
    dd= fasta_format(dd)
    selected=input$FastaList_rows_selected
    selected = selected[length(selected)]
    if(input$checkbox==TRUE){
      index = use_unique(dd,selected,indata)
    } else {
      index=grep(dd[(selected),"Accession"],indata$Accession,fixed=T)
    }
    indind=indata[index,]
    mm = ProteinPlotMat(dd[as.numeric(selected),"Sequence"],indind)
    info = summaryFunction(mm[[1]],indind,dd[(selected),"Accession"])
    return(info)
    },options = list(searching = FALSE,paging = FALSE))
  
  ##################################################################    
  ## modification summary Fasta List (tab 1)
  ##################################################################  

  
  mymod = reactive({
    dd= readIn2()
    indata = readIn1()
    indata=filter_amanda(indata,input$filter)
    dd= fasta_format(dd)
    selected=input$FastaList_rows_selected
    selected = selected[length(selected)]
})
  output$mod = renderDataTable({
    dd= readIn2()
    indata=readIn1()
    indata=filter_amanda(indata,input$filter)
    dd= fasta_format(dd)
    selected=input$FastaList_rows_selected
    selected = selected[length(selected)]
    if(input$checkbox==TRUE){
      index = use_unique(dd,selected,indata)
    } else {
      index=grep(dd[(selected),"Accession"],indata$Accession,fixed=T)
    }
    indind=indata[index,]
    mm = ProteinPlotMat(dd[as.numeric(selected),"Sequence"],indind)
    res=(modSummary(mm[[1]]))
    colnames(res) = c("Position","Modification","# Modifications","Procentage of all peptides")
    out <<- res
    outname <<- dd[(selected),"Accession"]
    return(res)
  })
 
  ##################################################################    
  ## info summary from Summary (tab 3)
  ##################################################################  
  
  output$info1 = renderDataTable({
    dd= readIn2()
    indata=readIn1()
    indata=filter_amanda(indata,input$filter)
    dd= fasta_format(dd)
    selected=input$summary_rows_selected
    selected = selected[length(selected)]
    selected=sub("\\|","\\\\|",selected)
    sel=(grep(selected,dd[,"Accession"]))
    if(input$checkbox==TRUE){
      index = use_unique(dd,sel,indata)
    } else {
      index=grep(dd[(sel),"Accession"],indata$Accession,fixed=T)
    }
    indind=indata[index,]
    mm = ProteinPlotMat(dd[(sel),"Sequence"],indind)
    info = summaryFunction(mm[[1]],indind,dd[(sel),"Accession"])
    return(info)
  },options = list(searching = FALSE,paging = FALSE))
  
  ##################################################################    
  ## modification summary Summary (tab 3)
  ##################################################################  
  
  
  output$mod1 = renderDataTable({
    dd= readIn2()
    indata=readIn1()
    indata=filter_amanda(indata,input$filter)
    dd= fasta_format(dd)
    selected=input$summary_rows_selected
    selected = selected[length(selected)]
    selected=sub("\\|","\\\\|",selected)
    sel=(grep(selected,dd[,"Accession"]))
    if(input$checkbox==TRUE){
      index = use_unique(dd,sel,indata)
    } else {
      index=grep(dd[(sel),"Accession"],indata$Accession,fixed=T)
    }
    indind=indata[index,]
    mm = ProteinPlotMat(dd[as.numeric(sel),"Sequence"],indind)
    res=(modSummary(mm[[1]]))
    colnames(res) = c("Position","Modification","# Modifications","Procentage of all peptides")
    out <<- res
    outname <<- selected
    return(res)
  })
   
##########################################
datasetInput <- out
    

output$downloadData1 <- downloadHandler(
  filename = function() { 
    paste(paste(outname,paste(input$filter,input$checkbox,sep="_"),sep="_"),'.csv', sep='') 
  },
  content = function(file) {
    write.csv(out, file)
  }
)

output$downloadData <- downloadHandler(
  filename = function() { 
    paste(paste(outname,paste(input$filter,input$checkbox,sep="_"),sep="_"),'.csv', sep='') 
  },
  content = function(file) {
    write.csv(out, file)
  }
)

outputOptions(output,'downloadData', suspendWhenHidden=FALSE)

})
