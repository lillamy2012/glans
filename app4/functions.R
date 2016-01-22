t2g = read.table("../data/blast/trans2gene.tab")
rownames(t2g)=t2g$V1
colnames(t2g)=c("transcript","gene","scaffold","start","end")

cmdCreate <- function(infile, outfile,type){
  if(type=="1"){
    return(paste("tblastx -db ../data/blast/merged.fa -query ", infile, " -outfmt 7",  " > ",
                 outfile, sep = ""))
  }
  if(type=="2")
    return(paste("tblastn -db ../data/blast/merged.fa -query ", infile, "-num_alignments 5 -outfmt 7",  " > ",
                 outfile, sep = ""))
}

createLink <- function(val) {
  sprintf('<a href=%s target="_blank" class="btn btn-primary">View</a>',val)
}

createFA <- function(id){
  outfile=paste("fasta",id,".fa",sep="")
  paste("pyfasta extract --header --fasta path2fasta", id, " > ",
        outfile, sep = "")
}
map2gene <- function(id){
  g = t2g[id,"gene"]
  g
}

mapTranscriptSeq <- function(id,outfile,extra){
  #print(id)
  gene= t2g[id,"gene"]
  g = t2g[id,c("scaffold","start","end")]
  g[,"start"] = g[,"start"]-extra
  g[,"end"] = g[,"end"]+extra
  comd = paste("samtools faidx ../data/blast/Marchantia_polymorpha.main_genome.scaffolds.fasta",paste(g[,1],paste(g[,2],g[,3],sep="-"),sep=":"))#,">",outfile) 
  return(list(comd,gene))
  }
