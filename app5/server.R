library(shiny)
library(reshape2)
source("functions.R")

shinyServer(function(input, output,session) {
  selected=NULL
#  use=NULL
#  keep=NULL
  ctouse=NULL
  rpkm=NULL
  coltoUse=NULL
  stList=NULL
  group=NULL
  stats=NULL
  #tot=NULL
  
  createSubData <- reactive({ ## remove outlier samples
    input$goButton
    outl = isolate(input$outlier)
    print("here")
    if(!is.null(outl)){
      torm = colnames(counts)%in%outl
      cols = cols[!torm,,drop=FALSE]
      counts = makeDataSet(counts,outl)
      genes = makeDataSet(genes,outl)
      counts = filterData(genes,counts,1)
    }
    ctouse <<- counts
    coltoUse <<- cols 
    rpkm <<- genes
  })
  
  doStats = reactive({
    createSubData()
    stList <<- runStats(ctouse,coltoUse)
  })
    
  
  makeGroupSet = function(){
    stIn = stList[[1]][which(stList[[2]][,"padj"]<input$padj),]
    groupL= groups(input$Set,input$Group,notinc,ofInt,stIn,input$minFC,input$minavFC)
    group <<- groupL[[1]]
    stats <<- groupL[[2]]
    }
  
  
  output$matplot = renderPlot({
    doStats()
    makeGroupSet()
    dt = rpkm[group,]
    #print(rownames(dt))
    selected <<- dt
    if(nrow(dt)==0)
      return(NULL)
    par(mar=c(10,5,1,1))
    if(input$transf)
      plotD = asinh(dt)
    if(!input$transf)
      plotD = dt
    toclick = melt(cbind(rownames(plotD),plotD))
    toclick$variable2 = as.numeric(toclick$variable)
    use<<-toclick
    matplot(t(plotD),type="l",lty=1,xaxt="n",yaxt="n",ylab="rpkm")
    axis(1,at=1:ncol(dt),colnames(dt),las=2)
    ii = seq.int(range(plotD)[1],range(plotD)[2],length.out =6)
    if(input$transf)
      lab= round(sinh(ii))
    if(!input$transf)
      lab = round(ii)
    axis(2,at=seq.int(range(plotD)[1],range(plotD)[2],length.out =6),labels = lab)
  })
  
  output$click_info <- renderPrint({
    nearPoints(use, input$plot_click, addDist = FALSE,xvar="variable2",yvar="value",threshold = 10)
  })
  
  output$tot_data = renderDataTable({
    if (nrow(selected)==0){
      return(NULL)
    }
    cbind(gene=rownames(selected),selected)
  })
  output$downloadData <- downloadHandler(
    filename = function() { 
      paste(paste("G",input$Group,"S",input$Set,"P",input$padj,"MFC",input$minavFC,"avFC",input$minFC,sep="_"),'csv', sep='.') 
    },
    content = function(file) {
      write.csv(selected, file)
    }
  )
  
  output$selc = renderPlot({
    makeGroupSet()
    
    
    #dt=groupData()
    #input$padj
    #input$Group
    #input$Set
    
    
    #input$minavFC
    #input$minFC
    #a = which(dt$one==1 & dt$mdist1>input$minFC & dt$avdist1>input$minavFC)
    #b = which(dt$two==1 & dt$mdist2>input$minFC & dt$avdist2>input$minavFC)
    #stat12 = length(intersect(a,b))
    #stat1o2 = length(union(a,b))
    #stat1 = length(setdiff(a,b))
    #stat2 = length(setdiff(b,a))
    #stats=c(stat1o2,stat12,stat1,stat2)
    
    bp=barplot(stats,col=c("lightblue","darkblue","lightyellow","lightgreen"),names.arg=c("at least one","in both","first only","second only"),ylab="nr genes")
    text(x=bp,y=max(stats)/2,labels =stats,cex=2)
  })
  
  
  
  
})