library(dplyr)
indata = read.csv("/Users/elin.axelsson/Desktop/MS_H31_H33_allmodifications.csv",skip=2,sep=";",dec=",")
tab_data = tbl_df(indata)

groups_l3 = distinct(select(tab_data,Accession,Sequence,Modifications))
groups_l2 = distinct(select(tab_data,Accession,Sequence))
groups_l1 = distinct(select(tab_data,Accession))

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
      perGr[[i]]=filter(tab_data,Accession==as.character(groups[i,1]) & Sequence==as.character(groups[i,2]) & Modifications==as.character(groups[i,3]))
    }
    if (level==2){
      perGr[[i]]=filter(tab_data,Accession==as.character(groups[i,1]) & Sequence==as.character(groups[i,2]))
    }
    if (level==1){
      perGr[[i]]=filter(tab_data,Accession==as.character(groups[i,1]))
    }
    
      groups[i,"numb"]=nrow(perGr[[i]])
    groups[i,paste("sample",1:4,sep="")]=apply(perGr[[i]][,4:7],2,function(x) sum(x=="X"))
  }
  groups$gr1 = rowSums(groups[,paste("sample",1:2,sep="")])
  groups$gr2 = rowSums(groups[,paste("sample",3:4,sep="")])
  groups$diff = abs(groups$gr1 -groups$gr2) 
  return(groups)
}





