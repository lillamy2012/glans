library(shiny)
library(reshape2)
source("functions.R")
#load("defaultStat.Rdata")

outputDir <- "/Users/elin.axelsson/glans/app5/responses"

### for the checkbox
saveData <- function(data) {
  data <- t(data)
  # Create a unique file name
  #fileName <- sprintf("%s_%s.csv", as.integer(Sys.time()), digest::digest(data))
  fileName <- "tomo.csv"
  # Write the file to the local system
  write.csv(
    x = data,
    file = file.path(outputDir, fileName), 
    row.names = FALSE, quote = TRUE
  )
}
removeData <-function(){
  tt= unlink(paste(outputDir,"tomo.csv",sep="/"))
}

loadData <- function() {
  # Read all the files into a list
  files <- list.files(outputDir, full.names = TRUE)
  data <- lapply(files, read.csv, stringsAsFactors = FALSE) 
  # Concatenate all data together into one data.frame
  data <- do.call(rbind, data)
  data
}
mych = loadData()

########################

shinyServer(function(input, output,session) {
 #set variable to NULL 
  selected=NULL
  use=NULL
  keep=NULL
  ctouse=NULL
  rpkm=NULL
  coltoUse=NULL
  stList=NULL
  group=NULL
  stats=NULL
  summlist=NULL
  
  ## remove outlier samples
  createSubData <- reactive({
    input$goButton # redo every time the ouliers been changed and go button clicked
    outl = isolate(input$outlier)
    print("here")
    if(!is.null(outl)){
      torm = colnames(counts)%in%outl
      cols = cols[!torm,,drop=FALSE]
      counts = makeDataSet(counts,outl)
      genes = makeDataSet(genes,outl)
      counts = filterData(genes,counts,1)
    }
    ctouse <<- counts # set variables use throughout
    coltoUse <<- cols 
    rpkm <<- genes
  })
  
  ## run DESeq2 to get statistics
  doStats = reactive({
    print("start")
    outl = isolate(input$outlier)
    ma = sapply(outls,function(x) identical(sort(x),sort(outl)))
    createSubData()
    input$goButton 
    if(length(outl)==0){
        print("using preCal")
        print("1")
        stList <<- ds[[1]]
    } else if(sum(ma)==1){
        print("using preCal")
        a = which(ma==TRUE)
        print(a)
        stList <<- ds[[a]]
      } else{
      stList <<- runStats(ctouse,coltoUse)
      }
    print("end")
  })
    
  ## use the group and set settings to get genes in group+set (plus numbers for barplot)
  calcGroupStat = reactive({
    print("once")
    input$goButton 
    outl = isolate(input$outlier)
    ma = sapply(outls,function(x) identical(sort(x),sort(outl)))
    
    if(length(outl)==0){
      print("using preCal")
      summList <<- cl[[1]]
    } else if(sum(ma)==1){
      print("using preCal")
      a = which(ma==TRUE)
      print(a)
      summList <<- cl[[a]]
    } else{
      stIn = stList[[1]][which(stList[[2]][,"padj"]<input$padj),]
      summList <<- calcStatsAll(notinc,ofInt,stIn)
    }
      print("once2")
  })
  
  CalcwitGr = reactive({
    print("per comb")
    input$goButton 
    groupL= groupsIm(input$Set,input$Group,input$minFC,input$minavFC,summList)
    group <<- groupL[[1]]
    stats <<- groupL[[2]]
  })
  
  rowSelect <- reactive({
    paste(sort(unique(input[["rows"]])),sep=',')
  })
  
  
  observe({
    #updateTextInput(session, "collection_txt", value = rowSelect() ,label = "Foo:" )
    if(length(rowSelect())>0){
      keep<<- 1
    }
    if(length(rowSelect())>0){ 
      nn = as.numeric(rowSelect())
      mych<<-rownames(selected)[nn]
      saveData(mych)
    } 
    if(length(rowSelect())==0 & !is.null(keep)){ # remove last 
      mych<<- NULL
      removeData()
    } 
    if(length(mych)>0){
      keep<<- 1
    }
    
  })
  
  
  
  output$matplot = renderPlot({
    doStats()
    calcGroupStat()
    CalcwitGr()
    #makeGroupSet()
    dt = rpkm[group,]
    selected <<- dt
    if(nrow(selected)==0)
      return(NULL)
    par(mar=c(10,5,1,1))
    if(input$transf)
      plotD = asinh(selected)
    if(!input$transf)
      plotD = selected
    toclick = melt(cbind(rownames(plotD),plotD))
    toclick$variable2 = as.numeric(toclick$variable)
    use<<-toclick
    matplot(t(plotD),type="l",lty=1,xaxt="n",yaxt="n",ylab="rpkm")
    axis(1,at=1:ncol(selected),colnames(selected),las=2)
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
    createSubData()
    doStats()
    calcGroupStat()
    CalcwitGr()
    #makeGroupSet()
    if (nrow(selected)==0){
      return(NULL)
    }
    addCheckboxButtons <- paste0('<input type="checkbox" name="row', rownames(selected),'" value="', rownames(selected), '">',"")
    precn = unique(unlist(mych))
    prec = which(rownames(selected)%in%precn)
    addCheckboxButtons[prec] <- paste0('<input type="checkbox" name="row', rownames(selected)[prec],'" value="', rownames(selected)[prec], '"checked=1"', '">',"")
    selected = cbind(gene=rownames(selected),selected)
    cbind(Pick=addCheckboxButtons, selected,id=1:nrow(selected))
  },escape = FALSE,options = list(bSortClasses = TRUE, aLengthMenu = c(5, 25, 50), iDisplayLength = 25),
  callback = "function(table) {
table.on('change.dt', 'tr td input:checkbox', function() {
setTimeout(function () {
Shiny.onInputChange('rows', $(this).add('tr td input:checkbox:checked').parent().siblings(':last-child').map(function() {
return $(this).text();
}).get())
}, 10); 
});
}")
  
  
  
  
  output$downloadData <- downloadHandler(
    filename = function() { 
      paste(paste("G",input$Group,"S",input$Set,"P",input$padj,"MFC",input$minavFC,"avFC",input$minFC,sep="_"),'csv', sep='.') 
    },
    content = function(file) {
      write.csv(selected, file)
    }
  )
  
  output$selc = renderPlot({
    #makeGroupSet()
    calcGroupStat()
    CalcwitGr()
    bp=barplot(stats,col=c("lightblue","darkblue","lightyellow","lightgreen"),names.arg=c("at least one","in both","first only","second only"),ylab="nr genes")
    text(x=bp,y=max(stats)/2,labels =stats,cex=2)
  })
  
  
  
  
})