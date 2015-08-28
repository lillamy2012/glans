library(shiny)
library(dplyr)

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
  print(class(groups))
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