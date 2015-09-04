library(shiny)
library(dplyr)

fasta_format = function(infile){
  #nrSamples = length(grep("Sample",colnames(infile)))
  nrSamples = (ncol(infile)-4)/2
  seqs = infile[,c("X",paste0("X.",1:(nrSamples-1)))]
  seqtouse=matrix(NA,nrow(seqs))
  for (i in 1:nrow(seqs)){
    if(length(setdiff(unique(c(seqs[i,])),"No coverage."))>1)
      print("pp")
    else{
      seqtouse[i] = unlist(setdiff(unique(c(seqs[i,])),"No coverage."))
  
    }
  } 
  results = data.frame(Accession=infile$Accession,Sequence=seqtouse)
  }



table_function = function(data,groups,level){
  perGr = list()
  groups$numb=NA
  groups$sample1=NA
  groups$sample2=NA
  groups$sample3=NA
  groups$sample4=NA
  groups$gr1=NA
  groups$gr2=NA
  for (i in 1:nrow(groups)){
    if (level==3){
      perGr[[i]]=filter(data,Accession==as.character(groups[i,1]) & Sequence==as.character(groups[i,2]) & Modifications==as.character(groups[i,3]))
    } else if (level==2){
      perGr[[i]]=filter(data,Accession==as.character(groups[i,1]) & Sequence==as.character(groups[i,2]))
    } else if (level==1){
      perGr[[i]]=filter(data,Accession==as.character(groups[i,1]))
    }
    
    groups[i,"numb"]=nrow(perGr[[i]])
    groups[i,paste("sample",1:4,sep="")]=apply(perGr[[i]][,4:7],2,function(x) sum(x=="X"))
  }
  
  groups$gr1 = rowSums(groups[,paste("sample",1:2,sep="")])
  groups$gr2 = rowSums(groups[,paste("sample",3:4,sep="")])
  groups$diff = abs(groups$gr1 -groups$gr2) 
  groups$p = apply(groups[,c("gr1","gr2")],1,function(x)  y = round(binom.test(ceiling(x[1]/2),ceiling((x[1]+x[2])/2))$p.value,5))
  groups=arrange(groups,desc(unq),p,desc(diff),desc(numb))
  
  #rem=c("numb","gr1","gr2")
  return(groups) #[,!colnames(groups)%in%rem])
}


stats_function = function(data,groups,level){
  perGr = vector("numeric",nrow(groups))
  for (i in 1:nrow(groups)){
    if (level==3){
      perGr[[i]]=nrow(filter(data,Accession==as.character(groups[i,1]) & Sequence==as.character(groups[i,2]) & Modifications==as.character(groups[i,3])))
    } else if (level==2){
      perGr[[i]]=nrow(filter(data,Accession==as.character(groups[i,1]) & Sequence==as.character(groups[i,2])))
    } else if (level==1){
      perGr[[i]]=nrow(filter(data,Accession==as.character(groups[i,1])))
    }
  }
  return(perGr)
}

filter_amanda=function(data,filt){
  data = filter(data,Amanda.Score>filt[1]) # & Amanda.Score<filt[2])
  return(data)
}


filter_probability=function(){
  print("ok")
}

group_function=function(data){  
  groups_l3 = distinct(select(data,Accession,Sequence,Modifications,unq))
  groups_l2 = distinct(select(data,Accession,Sequence,unq))
  groups_l1 = distinct(select(data,Accession,unq))
  return(list(groups_l1,groups_l2,groups_l3))
}


ord_list=function(data,id,levels){
  if(levels==1){
    toShow=filter(data,Accession ==id[["Acc"]])
  } else if (levels==2){
    toShow=filter(data,Accession ==id[["Acc"]],Sequence==id[["Seq"]])
  } else if (levels==3){
    toShow=filter(data,Accession ==id[["Acc"]],Sequence==id[["Seq"]],Modifications==id[["Mod"]])
  }
  return(toShow)
}
 
makeTab=function(indata){
  gr = sapply(strsplit(indata$Accession,";"),length)
  indata$unq = ifelse(gr>1,0,1)
  tab_data = tbl_df(indata)  
  return(tab_data)
} 

#indata = read.csv("/Users/elin.axelsson/Desktop/MS_H31_H33_allmodifications.csv",skip=2,sep=";",dec=",")
#gr = sapply(strsplit(indata$Accession,";"),length)
#indata$unq = ifelse(gr>1,0,1)
#tab_data = tbl_df(indata)

prop=function(indata){
sun = indata[,9:10]
sunsplit1 = strsplit(sun[,1],";")
un = unique(unlist(sunsplit1))
newmat = matrix(NA,nrow(indata),length(un))
colnames(newmat) = un
sunsplit2 = strsplit(sun[,2],";")
temp = lapply(sunsplit2,strsplit,":")


ids=lapply(temp,function(x) sapply(x,"[",1))
ps = lapply(temp,function(x) sapply(x,"[",2))
times = rep(1:length(temp),times=sapply(temp,length))
tot = data.frame(unlist(ids),as.numeric(unlist(ps)),times)
}

checkIfMap <- function(fasta,Sequences){
  match_ind=(str_locate_all(pattern = Sequences, fasta))
  wrong=which(is.na(sapply(match_ind,"[",1)))
  return(wrong)
}


ProteinPlotMat <- function(fasta,indind){
  by_character = strsplit(fasta,"")[[1]]
  covmat = matrix(0,length(by_character))
  Sequences = indind$Sequence
  wrongMap = checkIfMap(fasta,Sequences)
  if (length(wrongMap)==nrow(indind)){
    return(NULL)
  } else if(length(wrongMap)>0 ){
    indind=indind[-wrongMap,]
    Sequences = indind$Sequence
  } 
  Modifications = indind$Modifications
  Modifications=sub("N-Term",1,Modifications)
  match_ind=(str_locate_all(pattern = Sequences, fasta))
  modtype=regmatches(Modifications, gregexpr("(?<=\\().*?(?=\\))", Modifications, perl=T))
  print(modtype)
  for (i in 1:length(match_ind)){
    ii=match_ind[[i]]
    if(!is.na(ii[1]))
      covmat[ii[1]:ii[2]]=covmat[ii[1]:ii[2]]+1
  }
  
  position <- strsplit(Modifications, "[^[:digit:]]")
  v=list()
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
  cols =c("darkblue","maroon","pink","brown","yellow","orange","red","blue","green","mistyrose","lightblue","lightgrey")
  need= ncol(totmat)
  bp = barplot(t(totmat[,ncol(totmat):1]),col=cols[(length(cols)-need+1):length(cols)],border = NA,names.arg=rep("",nrow(totmat)))
  chars = rep("",nrow(totmat))
  isnot = which(rowSums(totmat[,2:ncol(totmat),drop=F])>0)
  chars[isnot]=by_character[isnot]
  mtext(side=1,at = bp,chars,cex=0.8)
  num=1:(length(by_character))
  num[!num%in%isnot]=""
  num = as.numeric(num)-1
  mtext(side=1,at = bp,num,cex=0.6,line=1,las=2)
  legend("topright",legend=colnames(totmat),fill=cols[length(cols):(length(cols)-need+1)])
  return(totmat)
}

summaryFunction <- function(mm,indind,accession){
  if(nrow(indind)==0){
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
        res[[k]] =(c(rownames(mm)[i],colnames(mm)[j],mm[i,j],(mm[i,j]/tot[i])*100))
    }
    }
    
  }
  re = do.call("rbind",res)
  return(re)
  }
  
tab=list(
"Oxidation"="no",
"Acetyl"="Ac",        
"Propionyl"="no",     
"Gln->pyro-Glu"="no",
"Ubi_LRGG"="UBi",     
"GlyGly"="UBi",      
"Trimethyl"="3meth",   
"Dimethyl"="2meth",    
"Methyl"= "Meth",       
"Phospho"= "Ph",    
"PropMeth"="Meth" ,  
"Methylthio"="no"
)
convert <-function(indata,table=tab){
for(i in 1:length(table)){
  #print(i)
  #print(names(table)[i])
  #print(table[[i]])
  #print(indata$Modifications)
  indata$Modifications=(gsub(names(table)[i],table[[i]],indata$Modifications))
}
  return(indata)
}  

#ttt = convert(indata,tab)  
#}
