library(shiny)
library(reshape2)

outputDir <- "/Users/elin.axelsson/glans/app5/responses"
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
print(paste(outputDir,"tomo.csv",sep="/"))
    print(tt)
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


shinyServer(function(input, output,session) {
  load("table.Rdata")
  tot = cbind(gene=rownames(tot),tot)
  selected=NULL
  use=NULL
  keep=NULL
  #updateTextInput(session, "collection_txt", value = tet ,label = "Foo:" )
  
  createSubSetData=function(){
    sel=tot[which(tot$group==input$Group),]
    if(input$Group==1){
      sel$one = sel$conditionMp_sperm_big
      sel$two = sel$conditionYoung_male_rep_big
      sel$avdist1 = sel$conditionMp_sperm_avdist
      sel$avdist2 = sel$conditionYoung_male_rep_avdist
      sel$mdist1 = sel$conditionMp_sperm_mdist
      sel$mdist2 = sel$conditionYoung_male_rep_mdist
    }
    if(input$Group==2){
      sel$one = sel$condition10DAF_big
      sel$two = sel$condition20DAF1mmembryo_big
      sel$avdist1 = sel$condition10DAF_avdist
      sel$avdist2 = sel$condition20DAF1mmembryo_avdist 
      sel$mdist1 = sel$condition10DAF_mdist
      sel$mdist2 = sel$condition20DAF1mmembryo_mdist
    }
    if(input$Group==3){
      sel$one = sel$v2_condition20DAF1mmembryo_big
      sel$two = sel$condition20DAFSoftembryo_big
      sel$avdist1 = sel$v2_condition20DAF1mmembryo_avdist
      sel$avdist2 = sel$condition20DAFSoftembryo_avdist
      sel$mdist1 = sel$v2_condition20DAF1mmembryo_mdist
      sel$mdist2 = sel$condition20DAFSoftembryo_mdist
    }
    
    if(input$Group==4){
      sel$one = sel$v2_condition20DAFSoftembryo_big
      sel$two = sel$conditionSpore_big
      sel$avdist1 = sel$v2_condition20DAFSoftembryo_avdist
      sel$avdist2 = sel$conditionSpore_avdist
      sel$mdist1 = sel$v2_condition20DAFSoftembryo_mdist
      sel$mdist2 = sel$conditionSpore_mdist
    }
    selected1 <<- sel
  }
  
  groupData = function(){
    sel = createSubSetData()
    sel = sel[sel$'res_anova$padj'<input$padj,]
    a = which(sel$one==1 & sel$mdist1>input$minFC & sel$avdist1>input$minavFC)
    b = which(sel$two==1 & sel$mdist2>input$minFC & sel$avdist2>input$minavFC)
    if(input$Set=="union")
      sel = sel[union(a,b),]
    if(input$Set=="intersect")
      sel = sel[intersect(a,b),]
    if(input$Set=="setdiff1")
      sel = sel[setdiff(a,b),]
    if(input$Set=="setdiff2")
      sel = sel[setdiff(b,a),]
    return(sel)
  }

  output$totdata=renderDataTable({
    dt=groupData()
    dt[,1:23]
    selected<<-dt[,1:23]
    dt
})
  rowSelect <- reactive({
    paste(sort(unique(input[["rows"]])),sep=',')
  })
  
  
  observe({
    updateTextInput(session, "collection_txt", value = rowSelect() ,label = "Foo:" )
    if(length(rowSelect())>0){
      keep<<- 1
    }
    if(length(rowSelect())>0){ 
      mych<<- as.numeric(rowSelect())
      saveData(mych)
    } 
    print(!is.null(keep))
    if(length(rowSelect())==0 & !is.null(keep)){ # remove last 
      mych<<- NULL
      removeData()
    } 
    if(length(mych)>0){
      keep<<- 1
    }
    
  })
  
  output$tot_data = renderDataTable({
    dt=groupData()
    dt[,1:23]
    selected<<-dt[,1:23]
    addCheckboxButtons <- paste0('<input type="checkbox" name="row', rownames(dt),'" value="', rownames(dt), '">',"")
    #Display table with checkbox buttons
    prec = unique(unlist(mych))
    addCheckboxButtons[prec] <- paste0('<input type="checkbox" name="row', rownames(dt)[prec],'" value="', rownames(dt)[prec], '"checked=1"', '">',"")
    cbind(Pick=addCheckboxButtons, dt,id=1:nrow(dt))
  },escape = FALSE, options = list(bSortClasses = TRUE, aLengthMenu = c(5, 25, 50), iDisplayLength = 25),
  callback = "function(table) {
    table.on('change.dt', 'tr td input:checkbox', function() {
  setTimeout(function () {
  Shiny.onInputChange('rows', $(this).add('tr td input:checkbox:checked').parent().siblings(':last-child').map(function() {
  return $(this).text();
  }).get())
  }, 10); 
});
  }")
  
 
  
  output$matplot = renderPlot({
    dt=groupData()
    dt =dt[,2:23]
    input$padj
    input$Type
    input$InType
    input$minavFC
    input$minFC
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
    axis(1,at=1:22,colnames(dt),las=2)
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
  
  output$selc = renderPlot({
    dt=groupData()
    input$padj
    input$Group
    input$Set
    
    
    input$minavFC
    input$minFC
    a = which(dt$one==1 & dt$mdist1>input$minFC & dt$avdist1>input$minavFC)
    b = which(dt$two==1 & dt$mdist2>input$minFC & dt$avdist2>input$minavFC)
    stat12 = length(intersect(a,b))
    stat1o2 = length(union(a,b))
    stat1 = length(setdiff(a,b))
    stat2 = length(setdiff(b,a))
    stats=c(stat1o2,stat12,stat1,stat2)
    bp=barplot(stats,col=c("lightblue","darkblue","lightyellow","lightgreen"),names.arg=c("at least one","in both","first only","second only"),ylab="nr genes")
    text(x=bp,y=max(stats)/2,labels =stats,cex=2)
    })
  
  output$downloadData <- downloadHandler(
    filename = function() { 
      paste(paste("G",input$Group,"S",input$Set,"P",input$padj,"MFC",input$minavFC,"avFC",input$minFC,sep="_"),'csv', sep='.') 
    },
    content = function(file) {
      write.csv(selected, file)
    }
  )
})