library(shiny)
library(DT)
source("functions.R")
hits=NULL
out=NULL
outname=NULL
shinyServer(function(input, output) {
  
    output$Results <- DT::renderDataTable({
      input$action
      progress <- shiny::Progress$new()
      progress$set(message = "Computing data", value = 0)
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progress$close())
      
      # Create a callback function to update progress.
      # Each time this is called:
      # - If `value` is NULL, it will move the progress bar 1/5 of the remaining
      #   distance. If non-NULL, it will set the progress to that value.
      # - It also accepts optional detail text.
      updateProgress <- function(value = NULL, detail = NULL) {
        if (is.null(value)) {
          value <- progress$getValue()
          value <- value + (progress$getMax() - value) / 5
        }
        progress$set(value = value, detail = detail)
      }
        data=(input$text)
        if (data!=""){
            con=file("temp.fa", open = "w")
            writeLines(c(">Q1",data),con=con)
            close(con)
            cmds <- cmdCreate("temp.fa","out.out",input$radio)
            system(cmds)
            my_table=as.data.frame(system2("bash", "../script/blast/grepCmd.sh",stdout=TRUE))
            my_table$link <- createLink(my_table[,1])
            ex = read.table("out.out")
            gene = map2gene(ex[,2])
            my_table=cbind(gene,ex[,-1])#,my_table$link)
            hits<<- my_table[,2]
            colnames(my_table)=c("gene_id","transcript_id","% identity","alignment length", "mismatches", "gap opens","q. start" ,"q. end","s. start" ,"s. end","evalue","bit score") #,"igv link")
        } else
            return(NULL)
        my_table
        },options = list(lengthMenu = c(15, 30, 50), pageLength = 15),escape=FALSE,
    callback = JS
        ('table.on("click.dt", "tr", function() {
        tabs = $(".tabbable .nav.nav-tabs li a");
        $(tabs[1]).click();})'))
    
    output$seq <- renderUI({
      sel = hits[as.numeric(input$Results_rows_selected[length(input$Results_rows_selected)])]
      seq = mapTranscriptSeq(sel,'seq.fasta',input$num)[[1]]
      gene = mapTranscriptSeq(sel,'seq.fasta',input$num)[[2]]
      ff = system(seq,intern = TRUE)
      headr = ff[1]
      fasta = ff[-1]
      outname <<- sub(">","",paste(headr = ff[1],gene,sep="_"))
      out <<- c(paste(headr,gene,sep="-"),fasta)
      HTML(c(out[1],paste(fasta,collapse = "<br>")))
      
    })
    output$downloadFasta <- downloadHandler(
      filename = function() { 
        paste(outname,"fasta",sep=".")
      },
      content = function(file) {
        write.table(out, file,quote = FALSE, row.names = FALSE, col.names = FALSE)
      }
    )
    
})
