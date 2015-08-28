library(shiny)
library(dplyr)
library(DT)
source("functions.R")
library(ggplot2)

###################################

#indata = read.csv("/Users/elin.axelsson/Desktop/MS_H31_H33_allmodifications.csv",skip=2,sep=";",dec=",")
#gr = sapply(strsplit(indata$Accession,";"),length)
#indata$unq = ifelse(gr>1,0,1)
#tab_data = tbl_df(indata)

shinyServer(function(input, output) {

  readIn <- reactive({
  inFile <- input$file1
  
  if (is.null(inFile))
    return(NULL)
  dd = read.csv(inFile$datapath, header=input$header, sep=input$sep, 
                quote=input$quote,skip=2,dec=",")
  return(dd)
  })
  
  
  output$contents <-renderTable({
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    dd = read.csv(inFile$datapath, header=input$header, sep=input$sep, 
             quote=input$quote,skip=2,dec=",")
    return(dd)
    })

  dataInput <- reactive({
    indata = readIn()
    gr = sapply(strsplit(indata$Accession,";"),length)
    indata$unq = ifelse(gr>1,0,1)
    tab_data = tbl_df(indata)
    tab_data = filter_amanda(tab_data,input$filter)
    ll = group_function(tab_data)
    groups_l3 = ll[[3]]
    groups_l2 = ll[[2]]
    groups_l1 = ll[[1]]
    if(input$radio==1){
      gro = groups_l1
    } else if(input$radio==2){
      gro = groups_l2
    } else if(input$radio==3){
      gro = groups_l3
    }
    tab = table_function(tab_data,gro,input$radio)
    rem=c("numb","gr1","gr2","unq")
    tab = tab[,!colnames(tab)%in%rem]
    return(tab)
  })
  
  tableInput <- reactive({
    indata = readIn()
    gr = sapply(strsplit(indata$Accession,";"),length)
    indata$unq = ifelse(gr>1,0,1)
    tab_data = tbl_df(indata)
    tab_data = filter_amanda(tab_data,input$filter)
    ll = group_function(tab_data)
    groups_l3 = ll[[3]]
    groups_l2 = ll[[2]]
    groups_l1 = ll[[1]]
    if(input$radio==1){
      gro = groups_l1
    } else if(input$radio==2){
      gro = groups_l2
    } else if(input$radio==3){
      gro = groups_l3
    }
    tab = stats_function(tab_data,gro,input$radio)
    return(tab)
  })
    
  output$hist = renderPlot({
   da = tableInput()
    da[da>100]=105
    hist(da,n=max(da),col="green")
  })
  
  
  output$test = DT::renderDataTable({
    dd= readIn()
    if (is.null(dd))
      return(NULL)
    gr = sapply(strsplit(as.character(dd$Accession),";"),length)
    dd$unq = ifelse(gr>1,0,1)
    ll = data.frame(Protein=unique(dd$Accession[which(dd$unq==1)]))
    ll
  })
  
  output$plot = renderPlot({
    dd = dataInput()
    dd= as.data.frame(dd)
    gr = sapply(strsplit(dd$Accession,";"),length)
    dd$unq = ifelse(gr>1,0,1)
    print(head(dd))
    hasZero = which(apply(dd[,c("sample1","sample3")],1,min)==0)
    sub = data.frame(sample1=c(dd$sample1),sample3=c(dd$sample3),unq=as.factor(c(dd$unq)))
    ggplot(sub,aes(sample1, sample3,color=unq)) + geom_point()+scale_x_log10()+scale_y_log10()+geom_abline()
  })
  
  output$click_info <- renderPrint({
    nearPoints(as.data.frame(dataInput()), input$plot_click, addDist = FALSE)
  })
  
  output$plot1 = renderPlot({
    dd = dataInput()
    dd= as.data.frame(dd)
    gr = sapply(strsplit(dd$Accession,";"),length)
    dd$unq = ifelse(gr>1,0,1)
    sub = data.frame(sample2=c(dd$sample2),sample4=c(dd$sample4),unq=as.factor(c(dd$unq)))
    ggplot(sub,aes(sample2, sample4,color=unq)) + geom_point()+scale_x_log10()+scale_y_log10()+geom_abline()
  })
  
  
  
  output$mytable = DT::renderDataTable({
    tab=dataInput()
    tab$Accession = paste0("<a href='#filtered_data'>", tab$Accession, "</a>")
    or = colnames(tab)
    print(tab)
    return(tab)
  }, options = list(lengthMenu = c(5, 30, 50), pageLength = 15),escape=FALSE,callback = JS(
    'table.on("click.dt", "tr", function() {
    tabs = $(".tabbable .nav.nav-tabs li a");
    $(tabs[4]).click();})'))
  
  

  output$filtered_data = renderDataTable({
    indata = readIn()
    gr = sapply(strsplit(indata$Accession,";"),length)
    indata$unq = ifelse(gr>1,0,1)
    tab_data = tbl_df(indata)
    selected <- input$mytable_rows_selected
    if (is.null(selected)){
      tab = tab_data
    } else {
      selected = selected[length(selected)]
      print(selected)
      id = as.data.frame(dataInput()) 
      print(id[as.numeric(selected),])
      tab=ord_list(tab_data,list(Acc=id[selected,"Accession"],Seq=id[selected,"Sequence"],Mod=id[selected,"Modifications"]),input$radio)
    }
  },options = list(lengthMenu = c(5, 30, 50), paging = FALSE),escape=FALSE)

}
  
)
    