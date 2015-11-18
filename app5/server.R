library(shiny)
library(reshape2)

shinyServer(function(input, output) {
  load("table.Rdata")
  tot = cbind(gene=rownames(tot),tot)
  selected=NULL
  use=NULL
  
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

  output$tot_data=renderDataTable({
    dt=groupData()
    print(nrow(dt))
    dt[,1:23]
    selected<<-dt[,1:23]
    dt
})
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
    print(str(toclick))
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
    print(stats)
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