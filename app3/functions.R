library(shiny)
library(dplyr)

##########################################
## functions for renaming and labeling 
##########################################

identifySumo = function(data){
  regmatches(data,(gregexpr("(?<=\\()GG[A-Z]+",data, perl=T)))="sumo"
  regmatches(data,(gregexpr("(?<=\\()[A-Z]+GG",data, perl=T)))="sumo"
  data
}

colorList = list(
  "no mod." = "lightgrey",
  "Ac" = "aquamarine",         
  "UBi" = "red",
  "Trimethyl" = "mistyrose",   
  "Dimethyl" = "pink",   
  "Meth" = "yellow",       
  "Ph" = "chocolate",
  "Cr" = "blue",
  "Ci" = "green",
  "Su" = "burlywood"
)

convert <-function(indata,table=tab){
  indata$Modifications=identifySumo(indata$Modifications)
  for(i in 1:length(table)){
    indata$Modifications=(gsub(names(table)[i],table[[i]],indata$Modifications))
    regmatches(indata$Modifications,(gregexpr(";.*(?<=\\(\\))", indata$Modifications, perl=T))) =""
    regmatches(indata$Modifications,(gregexpr(".*(?<=\\(\\))", indata$Modifications, perl=T))) =""
  }
  return(indata)
}  

tab=list(
  "Oxidation"="",
  "Acetyl"="Ac",        
  "Propionyl"="",     
  "Gln->pyro-Glu"="",
  "Ubi_LRGG"="UBi",     
  "GlyGly"="UBi",      
  #"Trimethyl"="3meth",   
  #"Dimethyl"="2meth",    
  "Methyl"= "Meth",       
  "Phospho"= "Ph",    
  "PropMeth"="Meth",   
  "Methylthio"="",
  "Crotonylation"="Cr",
  "Deamidated" = "Ci",
  "sumo" = "Su"
)

######################################################
## rearrange data frame formats
######################################################

fasta_format = function(infile){
  nrSamples = (ncol(infile)-4)/2
  seqs = infile[,c("X",paste0("X.",1:(nrSamples-1)))]
  seqtouse=matrix(NA,nrow(seqs))
  for (i in 1:nrow(seqs)){
    if(length(setdiff(unique(c(seqs[i,])),"No coverage."))>1)
      print("Something is wrong")
    else{
      seqtouse[i] = unlist(setdiff(unique(c(seqs[i,])),"No coverage."))
    }
  } 
  results = data.frame(Accession=infile$Accession,Sequence=seqtouse)
  }

#######################################################
## filter data, and check if map as assigned
#######################################################

filter_amanda=function(data,filt){
  data = filter(data,Amanda.Score>filt[1])
  return(data)
}

use_unique=function(data,index,indata){
  selected=sub("\\|","\\\\|",data[index,"Accession"])
  grepWord=paste("^","$",sep=selected)
  index = grep(grepWord,indata$Accession,fixed=F)
  return(index)
  }

checkIfMap <- function(fasta,Sequences){
  match_ind=(str_locate_all(pattern = Sequences, fasta))
  wrong=which(is.na(sapply(match_ind,"[",1)))
  return(wrong)
}

#########################################################
## make a matrix with coverage per position (one protein)
#########################################################

ProteinPlotMat <- function(fasta,indind){ ### fasta is the sequence of protein in question, indind is the data table matching 
  by_character = strsplit(fasta,"")[[1]] ### one column per character
  covmat = matrix(0,length(by_character)) ### matrix with 0 to start
  Sequences = indind$Sequence ### sequences reported to map protein
  wrongMap = checkIfMap(fasta,Sequences) ### check if correctly mapped
  if (length(wrongMap)==nrow(indind)){
    return(NULL)
  } else if(length(wrongMap)>0 ){
    indind=indind[-wrongMap,]
    Sequences = indind$Sequence
  } 
  Modifications = indind$Modifications ### modifications reported 
  Modifications=sub("N-Term",1,Modifications)  ### rename
  match_ind=(str_locate_all(pattern = Sequences, fasta)) ### start and end for each sequence to fasta file (list one item per sequence)
  for (i in 1:length(match_ind)){ ## for each sequence mapped to protein
    ii=match_ind[[i]] ## start and end in fasta
    if(!is.na(ii[1]))
      covmat[ii[1]:ii[2]]=covmat[ii[1]:ii[2]]+1 ## add one where cover
  }
  modtype=regmatches(Modifications, gregexpr("(?<=\\().*?(?=\\))", Modifications, perl=T))  ### type of modificantions per sequence (list)
  position <- strsplit(Modifications, "[^[:digit:]]") # position in peptid of modification (list matching modtype)
  v=list() ## list with new position (relative fasta file)
  for(i in 1:length(position)){
    v[[i]]=as.numeric(unlist(position[i]))+unlist(match_ind[i])[1]-1 # 
  }
  protein_pos=lapply(v,function(x) x[!is.na(x)])
  tot = list()
  for(i in 1:length(v)){
    tot[[i]]= data.frame(modtype[[i]],protein_pos[[i]])
  }
  totdf = do.call("rbind",tot)
  if(nrow(totdf)>0){
    stat_nr = as.matrix(table(totdf))
    totmat = matrix(0,nrow(covmat),nrow(stat_nr)+1)
    totmat[,1]=covmat
    for(i in 1:nrow(stat_nr)){
      totmat[as.numeric(colnames(stat_nr)),i+1]=stat_nr[i,]
    }
    totmat[,1]=totmat[,1]-rowSums(totmat[,2:ncol(totmat),drop=F])
    colnames(totmat)=c("no mod.",rownames(stat_nr))
    rownames(totmat)=paste0(by_character,0:(length(by_character)-1))
  } else {
    totmat = covmat
  }
  return(list(totmat,by_character))
}

ProteinPlot <- function(totmat,by_character){
  cols =c("darkblue","darkred","pink","darkorchid","orange","red","blue","green","mistyrose","lightblue","lightgrey")
  cols2= c(unlist(colorList[colnames(totmat)]))
  extra2 = setdiff(colnames(totmat),names(colorList))
  if(length(extra2)>0){
    extra = setdiff(cols,cols2)
    colsTot = c(cols2,extra[1:length(extra2)])
    names(colsTot) = c(names(cols2),extra2)
  } else 
    colsTot=cols2
  bp = barplot(t(totmat[,c(names(colsTot)[c(2:length(colsTot),1)])]),col=colsTot[c(2:length(colsTot),1)],border = NA,names.arg=rep("",nrow(totmat)))
  chars = rep("",nrow(totmat))
  isnot = which(rowSums(totmat[,2:ncol(totmat),drop=F])>0)
  chars[isnot]=by_character[isnot]
  mtext(side=1,at = bp,chars,cex=0.8)
  num=1:(length(by_character))
  num[!num%in%isnot]=""
  num = as.numeric(num)-1
  mtext(side=1,at = bp,num,cex=0.6,line=1,las=2)
  legend("topright",legend=names(colsTot),fill=colsTot)
  return(totmat)
}

#####################################################################
## summarize
#####################################################################

summaryFunction <- function(mm,indind,accession){
  if(is.null(mm)){
    return(NULL)
  } else if(nrow(indind)==0){
    mysum = data.frame(coverage=0,nr.peptides=0,modification_types=0,modifications_detected=0)
   } else {
    cov=sum(mm[,1]>0)/sum(mm[,1]>-1)
    counts=nrow(indind)
    modif=ncol(mm)-1
    modif.counts=sum(mm[,-1])
    mysum = data.frame(coverage=cov,nr.peptides=counts,modification_types=modif,modifications_detected=modif.counts)
    rownames(mysum)=accession
   }
  return(mysum)  
}

modSummary <- function(mm){
  mod = rowSums(mm[,-1,drop=F])
  tot = rowSums(mm)
  res = list()
  k = 0 
  for (i in 1:nrow(mm)){
    for (j in 2:ncol(mm)){
      if(mm[i,j]>0){
        k = k+1
        res[[k]] =(c(rownames(mm)[i],colnames(mm)[j],mm[i,j],round((mm[i,j]/tot[i])*100,1)))
        }
      }
    }
  re = do.call("rbind",res)
  return(re)
  }

getSampleName <- function(infile){
  nrSamples = (ncol(infile)-4)/2
  ids = c("X$",paste0("X.",1:(nrSamples-1)))
  ind = sapply(ids,grep,colnames(infile))-1
  return(colnames(infile)[ind])
}
       
splitToGroups <- function(infile,group1,group2){
  #print("so far")
  #print(group1)
  gr1 = infile[which(infile[,group1]=="X"),]
  #print(head(gr1))
  gr2 = infile[which(infile[,group2]=="X"),]
  res = list(gr1,gr2)
  return(res)
}
 
  
