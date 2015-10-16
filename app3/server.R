options(stringsAsFactors=FALSE)
library(shiny)
library(dplyr)
library(DT)
source("functions.R")
library(stringr)
###################################
out=NULL
infile1=NULL
#infile1.1=NULL
#infile1.2=NULL
infile2=NULL
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
    infile1 <<- dd
    return(dd)
  })
  
  readIn2 <- reactive({
    inFile <- input$file2
    if (is.null(inFile))
      return(NULL)
    dd = read.csv(inFile$datapath, header=TRUE, sep=";", 
                  quote="'",skip=2,dec=",",comment="",stringsAsFactors=FALSE)
    dd = dd[-1,]
    dd = dd[!is.na(dd[,1]),]
    infile2 <<- dd
    return(dd)
  })
  
  
  ##################################################################    
  ## summarize proteins in study (tab 'Fasta List'). This is the first thing that happens, infile1 and infile2 are created
  ##################################################################  
  
  output$FastaList = DT::renderDataTable({
    if (is.null(infile1))
      infile1= readIn1()
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
## summarize peptides, modifications etc per protein (tab 'Summary')
##################################################################  
  
  output$summary =  DT::renderDataTable({
    dd = infile2
    grp=0
    indata=infile1
    res = matrix(NA,nrow=nrow(dd),ncol=4)
    dd= fasta_format(dd)
     if(is.null(indata)|is.null(dd))
      return(NULL)
    indata=filter_amanda(indata,input$filter)
    if(!is.null(input$group1) & !is.null(input$group2)){
      groups=splitToGroups(indata,input$group1,input$group2)
      infile1.1 = groups[[1]]
      infile1.2 = groups[[2]]
      grp=1
      res1 = matrix(NA,nrow=nrow(dd),ncol=4)
      res2 = res1
      res = list(res1,res2)
    }
    withProgress(message = 'Calculating', value = 0, {
      for (i in 1:nrow(dd)){
        fasta=dd$Sequence[i]
        if(input$checkbox==TRUE){
          selected=sub("\\|","\\\\|",dd[i,"Accession"])
          grepWord=paste("^","$",sep=selected)
          if(grp==1){
            index1 = grep(grepWord,infile1.1$Accession,fixed=F,value=F)
            index2 = grep(grepWord,infile1.2$Accession,fixed=F,value=F)
            indind1 = infile1.1[index1,]
            indind2 = infile1.2[index2,]
          } else {
            index1 = grep(grepWord,indata$Accession,fixed=F,value=F)
            indind1=indata[index1,]
          }
        } else {
          if(grp==1){
            index1 = grep(dd[i,"Accession"],infile1.1$Accession,fixed=F,value=F)
            index2 = grep(dd[i,"Accession"],infile1.2$Accession,fixed=F,value=F)
            indind1 = infile1.1[index1,]
            indind2 = infile1.2[index2,]
            
          } else {
            index1=grep(dd[i,"Accession"],indata$Accession,fixed=T)
            indind1=indata[index1,]
          }}
        
        if(class(res)=="list"){
          if(length(index1)>0){  
            res[[1]][i,]=unlist(summaryFunction(ProteinPlotMat(fasta,indind1)[[1]],indind1,dd[i,"Accession"]))
          } else 
            res[[1]][i,]=rep(0,ncol(res[[1]]))
          if(length(index2)>0){  
            res[[2]][i,]=unlist(summaryFunction(ProteinPlotMat(fasta,indind2)[[1]],indind2,dd[i,"Accession"]))
          } else 
            res[[2]][i,]=rep(0,ncol(res[[2]]))
         } else
          if(length(index1)>0){  
            res[i,]=unlist(summaryFunction(ProteinPlotMat(fasta,indind1)[[1]],indind1,dd[i,"Accession"]))
          } else {
            res[i,]=rep(0,ncol(res))
          }
      # Increment the progress bar, and update the detail text.
        incProgress(1/nrow(dd), detail = paste("Protein", dd[i,"Accession"]))
      # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
      }
    })
    if(class(res)=="list")
      res = cbind(res[[1]],res[[2]]) 
    res = as.data.frame(res)
    rownames(res)=dd$Accession
    if (grp==1)
      colnames(res)=c("coverage_gr1","nr.peptides_gr1","modification_types_gr1","modifications_detected_gr1",
                    "coverage_gr2","nr.peptides_gr2","modification_types_gr2","modifications_detected_gr2")
    if (grp==0)
    colnames(res)=c("coverage","nr.peptides","modification_types","modifications_detected")
    res
  },options = list(lengthMenu = c(15, 30, 50), pageLength = 15),escape=FALSE,callback = JS(
    'table.on("click.dt", "tr", function() {
    tabs = $(".tabbable .nav.nav-tabs li a");
    $(tabs[3]).click();})'))

  
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
    grp=0
    indata=readIn1()
    indata=filter_amanda(indata,input$filter)
    dd= fasta_format(dd)
    if(!is.null(input$group1) & !is.null(input$group2)){
      groups=splitToGroups(indata,input$group1,input$group2)
      infile1.1 = groups[[1]]
      infile1.2 = groups[[2]]
      grp=1
      res = list(infile1.1,infile1.2)
    } else 
      res = indata
    selected=input$summary_rows_selected
    selected = selected[length(selected)]
    selected=sub("\\|","\\\\|",selected)
    sel=(grep(selected,dd[,"Accession"]))
    if(input$checkbox==TRUE){
      if(grp==1){
        index1 = use_unique(dd,sel,res[[1]])
        index2 = use_unique(dd,sel,res[[2]])
        indind1=indata[index1,]
        indind2=indata[index2,]
      }else  {
        index = use_unique(dd,sel,res)
        indind=indata[index,]
      }
    }
    if(input$checkbox==FALSE)
      if(grp==1){
        index1=grep(dd[(sel),"Accession"],res[[1]]$Accession,fixed=T)
        index2=grep(dd[(sel),"Accession"],res[[2]]$Accession,fixed=T)
        indind1=res[[1]][index1,]
        indind2=res[[2]][index2,]
    } else {
      index=grep(dd[(sel),"Accession"],res$Accession,fixed=T)
      indind=indata[index,]
    }
   if(grp==1){
     if(nrow(indind1)>0){
      mm1 = ProteinPlotMat(dd[as.numeric(sel),"Sequence"],indind1)
     }
     if(nrow(indind2)>0){
       mm2 = ProteinPlotMat(dd[as.numeric(sel),"Sequence"],indind2)
     } 
     par(mfrow=c(2,1))
     ProteinPlot(mm1[[1]],mm1[[2]])
     ProteinPlot(mm2[[1]],mm2[[2]])
     
   } else {
      if(nrow(indind)>0){
        mm = ProteinPlotMat(dd[as.numeric(sel),"Sequence"],indind)
        ProteinPlot(mm[[1]],mm[[2]])
    }
    else
      plot(1:100,pch=NULL)
  }})
    
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
    colnames(res) = c("Position","Modification","# Modifications","Percentage of all peptides")
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
    colnames(res) = c("Position","Modification","# Modifications","Percentage of all peptides")
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

output$done <- reactive({
  dd= readIn2()
  samp = getSampleName(dd)
  oo <<-samp
  colnames(dd)
  #print(dd)
  indata=filter_amanda(infile1,input$filter)
  #ProtPerSample(indata,dd)
})

output$group <- reactive({
  if(input$grouping==FALSE){
    dd= readIn2()
    samp = getSampleName(dd)
    oo <<-samp
    yes=TRUE
  }
  else{
    yes=FALSE
  oo=NULL
  }
    })
  
outputOptions(output,'done', suspendWhenHidden=FALSE)
outputOptions(output,'group', suspendWhenHidden=FALSE)

output$oo <-renderText({
  dd= readIn2()
  samp = getSampleName(dd)
  oo=samp
})

output$group1 <- renderUI({
  dd= readIn2()
  samp = getSampleName(dd)
  checkboxGroupInput("group1", "Group1", 
                     choices  = samp,
                     selected = NULL)
})
output$group2 <- renderUI({
  dd= readIn2()
  samp = getSampleName(dd)
  checkboxGroupInput("group2", "Group2", 
                     choices  = samp,
                     selected = NULL)
})


output$plotGroups <- renderPlot({
  gr1 = input$group1
  print(head(gr1))
  gr2 = input$group2
  indata=filter_amanda(infile1,input$filter)
  data= ProtPerSample(indata,infile2)
  newdata=data/colSums(data)
  plot(rowSums(newdata[,gr1,drop=F]),rowSums(newdata[,gr2,drop=F]),ylab="group2",xlab="group1")
  abline(0,1,col="red",lwd=3)
  
})
})