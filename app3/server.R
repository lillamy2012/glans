options(stringsAsFactors=FALSE)
library(shiny)
library("Hmisc")
library(dplyr)
library(DT)
source("functions.R")
library(stringr)

###################################
ss=NULL
leS=0
leF=0
out=NULL
out2=NULL
outname=NULL
info=NULL
res=NULL
res2=NULL
infile1=NULL
infile2=NULL
toclick=NULL
summ=NULL
returnpdf=FALSE
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
    dd
  },options = list(lengthMenu = c(5, 30, 50), pageLength = 15),escape=FALSE,callback = JS(
    'table.on("click.dt", "tr", function() {
    tabs = $(".tabbable .nav.nav-tabs li a");
    $(tabs[3]).click();})'))
  
  ##################################################################    
  ## plot protein from fasta (tab 1)
  ##################################################################  
  
  
  
  
  
##################################################################    
## summarize peptides, modifications etc per protein (tab 'Summary')
##################################################################  
  
  output$summary =  DT::renderDataTable({
    summ=NULL
    indata=infile1 
    dd = infile2
    if(is.null(indata)|is.null(dd))
      return(NULL)
    dd= fasta_format(dd)
    indata=filter_amanda(indata,input$filter) ## amanda score
    if(input$checkbox==TRUE) ## unique only
      indata = indata[grep(";",indata$Accession,invert=TRUE),]
    data= ProtPerSample(indata,infile2)
    input$goButton
    isolate({
    g1 = input$group1
    g2 = input$group2
    if(!is.null(g1) & !is.null(g2)){ ## pairwise
      groups=isolate(splitToGroups(indata,input$group1,input$group2))
      infile1.1 = groups[[1]]
      infile1.2 = groups[[2]]
      grp=1
      tot1 = sum(data[,input$group1])
      tot2 = sum(data[,input$group2])
      res1 = matrix(NA,nrow=nrow(dd),ncol=5)
      res2 = res1
      res = list(res1,res2)
    } else { # not pairwise
      grp=0
      res = matrix(NA,nrow=nrow(dd),ncol=5)
      tot1 = sum(data)
    }
    })
    withProgress(message = 'Calculating', value = 0, {
      for (i in 1:nrow(dd)){
        fasta=dd$Sequence[i]
        grepWord = dd[i,"Accession"]
          if(grp==1){ # pairwise
            index1 = grep(grepWord,infile1.1$Accession,fixed=F,value=F)
            index2 = grep(grepWord,infile1.2$Accession,fixed=F,value=F)
            indind1 = infile1.1[index1,]
            indind2 = infile1.2[index2,]
            if(length(index1)>0 & !is.null(ProteinPlotMat(fasta,indind1))){  
              res[[1]][i,]=c(unlist(summaryFunction(ProteinPlotMat(fasta,indind1)[[1]],indind1,dd[i,"Accession"])),tot1)
            } else 
              res[[1]][i,]=c(rep(0,ncol(res[[1]])-1),tot1)
            if(length(index2)>0 & !is.null(ProteinPlotMat(fasta,indind2))){  
              res[[2]][i,]=c(unlist(summaryFunction(ProteinPlotMat(fasta,indind2)[[1]],indind2,dd[i,"Accession"])),tot2)
            } else 
              res[[2]][i,]=c(rep(0,ncol(res[[2]])-1),tot2)
          }
        
          if (grp==0) { # not pairwise
            index1 = grep(grepWord,indata$Accession,fixed=F,value=F)
            indind1=indata[index1,]
            if(length(index1)>0 & !is.null(ProteinPlotMat(fasta,indind1))){  
              res[i,]=c(unlist(summaryFunction(ProteinPlotMat(fasta,indind1)[[1]],indind1,dd[i,"Accession"])),tot1)
            } else {
              res[i,]=c(rep(0,ncol(res)-1),tot1)
            }
          }
        
      # Increment the progress bar, and update the detail text.
        incProgress(1/nrow(dd), detail = paste("Protein", dd[i,"Accession"]))
      # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
      }
    })
    if(grp==1){
      res = cbind(res[[1]],res[[2]]) 
    }
    res = as.data.frame(res)
    rownames(res)=dd$Accession
    if (grp==1)
      colnames(res)=c("coverage_gr1","nr.peptides_gr1","modification_types_gr1","modifications_detected_gr1","tot1",
                    "coverage_gr2","nr.peptides_gr2","modification_types_gr2","modifications_detected_gr2","tot2")
    if (grp==0)
      colnames(res)=c("coverage","nr.peptides","modification_types","modifications_detected","tot")
    summ <<- res
    res
  },options = list(lengthMenu = c(15, 30, 50), pageLength = 15),escape=FALSE,callback = JS(
    'table.on("click.dt", "tr", function() {
    tabs = $(".tabbable .nav.nav-tabs li a");
    $(tabs[3]).click();})'))

  
  
  ##################################################################    
  ## plot protein from Summary (tab 3)
  ##################################################################  
  
  output$my_prot1 =renderPlot({
    info <<- NULL
    res <<- NULL
    res2 <<-NULL
    dd= readIn2()
    indata=readIn1()
    if(is.null(indata)|is.null(dd))
      return(NULL)
    dd= fasta_format(dd)
    track1 = track()
    sel=ss
    if(length(sel)<1)
      return(NULL)
    indata=filter_amanda(indata,input$filter) # filter on amanda score
    if(input$checkbox==TRUE){ # only want to use unique peptides 
      multi=grep(";",indata$Accession)
      uniq= grep(";",indata$Accession,invert=TRUE)
      indata=indata[uniq,]
      }
    data= ProtPerSample(indata,infile2)
    if(!is.null(input$group1) & !is.null(input$group2)){ # paired analysis
      groups=(splitToGroups(indata,input$group1,input$group2))
      grp=1
      res=list()
      tot=c(1,1)
      res[[1]] = groups[[1]]
      res[[2]] = groups[[2]]
      if(input$norm){
        tot[1]=sum(data[,input$group1])
        tot[2]=sum(data[,input$group2])
      }
    } else { # combine all
      res = indata
      if(input$norm){
        tot = sum(data)
      } else 
        tot=1
      grp=0
    }
    if(grp==1){
      index1=grep(dd[(sel),"Accession"],res[[1]]$Accession,fixed=T)
      index2=grep(dd[(sel),"Accession"],res[[2]]$Accession,fixed=T)
      indind1=res[[1]][index1,]
      indind2=res[[2]][index2,]
      if(nrow(indind1)>0){
        mm1 = ProteinPlotMat(dd[as.numeric(sel),"Sequence"],indind1)
        info1 = summaryFunction(mm1[[1]],indind1,dd[(sel),"Accession"])
        res <<- modSummary(mm1[[1]])
      } else{
          mm1 = list(matrix(0,nrow=1,ncol=1),NULL)
      }
      if(nrow(indind2)>0){
        mm2 = ProteinPlotMat(dd[as.numeric(sel),"Sequence"],indind2)
        info2 = summaryFunction(mm2[[1]],indind2,dd[(sel),"Accession"])
        res2 <<- modSummary(mm2[[1]])
      } 
      info_t = rbind(info1,info2)
      info_t$id = c(rownames(info1),rownames(info2))
      info_t = info_t[,c(5,1:4)]
      rownames(info_t)=c("group1","group2")
      info <<- info_t
      par(mfrow=c(2,1))
      ys = max(max(mm2[[1]]/tot[2]),max(mm1[[1]]/tot[1]))
      ProteinPlot(mm1[[1]]/tot[1],mm1[[2]],ylim=c(0,ys),input$returnpdf,input$size,type=2)
      ProteinPlot(mm2[[1]]/tot[2],mm2[[2]],ylim=c(0,ys),input$returnpdf,input$size,type=3)
        
    }
    if(grp==0) {
      index=grep(dd[(sel),"Accession"],res$Accession,fixed=T)
      indind=indata[index,]
      if(nrow(indind)>0){
        mm = ProteinPlotMat(dd[as.numeric(sel),"Sequence"],indind)
        info <<- summaryFunction(mm[[1]],indind,dd[(sel),"Accession"])
        res <<- modSummary(mm[[1]])
        ProteinPlot(mm[[1]]/tot,mm[[2]],input$returnpdf,input$size,type=1)
        
      }
    }
  })
  ##################################################################    
  ## modification summary Summary (tab 3)
  ##################################################################  
  
  output$mod1 = renderDataTable({
    track1 = track()
    sel=ss
    input$summary_rows_selected
    input$FastaList_rows_selected
    input$filter
    input$group1
    input$group2
    input$checkbox
    if(is.null(res))
      return(NULL)
    if(length(sel)<1)
      return(NULL)
    colnames(res) = c("Position","Modification","# Modifications","Percentage of all peptides")
    out <<- res
    outname <<- infile2[sel,"Accession"]
    return(res)
  })
  
  output$mod2 = renderDataTable({
    track1 = track()
    sel=ss
    input$summary_rows_selected
    input$FastaList_rows_selected
    input$filter
    input$group1
    input$group2
    input$checkbox
    if(is.null(res2)){
      return(NULL)
    }
    if(length(sel)<1)
      return(NULL)
    colnames(res2) = c("Position","Modification","# Modifications","Percentage of all peptides")
    out2 <<- res2
    outname <<- infile2[sel,"Accession"]
    return(res2)
  })
  
  ########################################################
  ## download tab 3
  #########################################################
  output$downloadData0 <- downloadHandler(
    filename <- function(){
      paste(paste(outname,paste("a",input$filter,"u",input$checkbox,"n",input$norm,"g",!input$grouping,input$group1,sep="_"),sep="_"),'.pdf', sep='')
    },
    content = function(file) {
     # if(!is.null(input$group1)){
      file.copy("plot.pdf", file)
    #  } else {
     # file.copy("tmp.pdf", file)
      #}
    }
  )
  
  output$downloadData1 <- downloadHandler(
    filename = function() { 
      paste(paste(outname,paste("a",input$filter,"u",input$checkbox,"n",input$norm,"g",!input$grouping,input$group1,sep="_"),sep="_"),'.csv', sep='') 
    },
    content = function(file) {
      write.csv(out, file)
    }
  )
  
  output$downloadData2 <- downloadHandler(
    filename = function() { 
      paste(paste(outname,paste("a",input$filter,"u",input$checkbox,"n",input$norm,"g",!input$grouping,input$group2,sep="_"),sep="_"),'.csv', sep='') 
    },
    content = function(file) {
      write.csv(out2, file)
    }
  )
  
 output$downloadSummary <- downloadHandler(
   filename = function() { 
     paste(paste("a",input$filter,"u",input$checkbox,"n",input$norm,"g",!input$grouping,input$group1,input$group2,sep="_"),'.csv',sep="")
   },
   content = function(file) {
     write.csv(summ, file)
   }
 )
  
##################################################################    
## info summary from Summary (tab 3)
##################################################################  
  
  output$info1 = renderDataTable({
    input$summary_rows_selected
    input$FastaList_rows_selected
    input$filter
    input$group1
    input$group2
    input$checkbox
    return(info)
  },options = list(searching = FALSE,paging = FALSE))
  
  ##################################################################    
  ## Scatter plot (tab 5)
  ##################################################################  
  
  output$plotGroups <- renderPlot({
    gr1 = input$group1
    gr2 = input$group2
    input$norm
    if(is.null(gr1) | is.null(gr2)){
      plot(1,1,col="white",axes=FALSE,ylab="",xlab="")
      text(1,1,"Need to define to groups before you can see this plot",cex=1)
    } else{
      indata=filter_amanda(infile1,input$filter)
      if(input$checkbox==TRUE) ## unique only
        indata = indata[grep(";",indata$Accession,invert=TRUE),]
      if(input$norm){
        data= ProtPerSample(indata,infile2)
        pp1 = getConfInt(data[,gr1,drop=F])
        pp2 = getConfInt(data[,gr2,drop=F])
        df = data.frame(gr1=pp1[1,],gr2=pp2[1,])
        toclick <<-df 
        plot(x=pp1[1,],y=pp2[1,],asp=1,pch=19,ylab="group2",xlab="group1")
        segments(pp1[2,],pp2[1,],pp1[3,],pp2[1,])
        segments(pp1[1,],pp2[2,],pp1[1,],pp2[3,])
        abline(0,1,col="red",lwd=3)
      } else{
        plot(1,1,col="white",axes=FALSE,ylab="",xlab="")
        text(1,1,"This type of plot doesn't make sense when data not is normalized",cex=1)
      }
    } 
    
  })
  output$click_info <- renderPrint({
    nearPoints(toclick, input$plot_click, addDist = FALSE,xvar="gr1",yvar="gr2")
  })
  


##################################################
## help functions
##################################################
  
## has file2 been uploaded?
output$done <- reactive({
  dd= readIn2()
  if(is.null(dd))
    return(NULL)
  return("yes")
})

## grouped analysis or not?
output$group <- reactive({
  if(input$grouping==FALSE){
    dd= readIn2()
    yes=TRUE
  }
  else{
    yes=FALSE
  }
    })
  
## to show boxes always
outputOptions(output,'done', suspendWhenHidden=FALSE)
outputOptions(output,'group', suspendWhenHidden=FALSE)


## choose group1
output$group1 <- renderUI({
  dd= readIn2()
  samp = getSampleName(dd)
  checkboxGroupInput("group1", "Group1", 
                     choices  = samp,
                     selected = NULL)
})
## choose group2
output$group2 <- renderUI({
  dd= readIn2()
  samp = getSampleName(dd)
  checkboxGroupInput("group2", "Group2", 
                     choices  = samp,
                     selected = NULL)
})

track <- function(){
  ## no selctions at all
  if(is.null(input$FastaList_rows_selected) & is.null(input$summary_rows_selected)){
    new=NULL
    print("both null") 
    ss <-- new
    return(NULL)
  }
  dd=infile2
  if(is.null(dd))
    return(NULL)
  ## set up fl
  if(!is.null(input$FastaList_rows_selected)){
    fl1 = as.numeric(input$FastaList_rows_selected[length(input$FastaList_rows_selected)])
  }
  ## set up summary
  if(!is.null(input$summary_rows_selected)){
    selected=input$summary_rows_selected
    selected = selected[length(selected)]
    selected=sub("\\|","\\\\|",selected) # id
    su1=(grep(selected,dd[,"Accession"])) # 
  }
  leFasta = length(input$FastaList_rows_selected)
  leSummary = length(input$summary_rows_selected)
  if(leSummary<leS | leFasta<leF){
    print("disselection")
    leS<<-leSummary
    leF<<-leFasta
    new = ss
  }
  if(leSummary>leS){
    print("summary selction")
    leS <<-leSummary
    new = su1
    ss <<- new
  }
  
  if(leFasta>leF){
    print("fasta selection")
    leF <<-leFasta
    new = fl1
    ss <<- new
  }
  return(new)
}
  })
