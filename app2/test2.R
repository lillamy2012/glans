library(stringr)
indata = read.csv("/Users/elin.axelsson/Desktop/MS_H31_H33_allmodifications.csv",skip=2,sep=";",dec=",")
fasta = "MARTKQTARKSHGGKAPRTLLATKAARKSAPTTGGVKKPHRYRPGTVALREIRKYQKSTELLIRKLPFQRLVREIAQDYKTDLRFQSHAVLALQEAAEAYLVGLFEDTNLCAIHAKRVTIMPKDVQLARRIRGERA*"

#u1 = "AARKSAPTTGGV"
#indind = indata[which(indata$Accession=="AT3G27360|H3.1"),]
#indind = indata[which(indata$Accession=="AT4G27230|H2A.2"),]
index=grep("AT3G27360\\|H3.1",indata$Accession)
indind=indata[index,]
fasta="MARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRFRPGTVALREIRKYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSSAVAALQEAAEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA*"
#test = "MAGRGKQLGSGAAKKSTSRSSKAGLQFPVGRIARFLKAGKYAERVGAGAPVYLAAVLEYLAAEVLELAGNAARDNKKTRIVPRHIQLAVRNDEELSKLLGDVTIANGGVMPNIHNLLLPKKAGSSKPTEED*"


ProteinPlot <- function(fasta,indind){
  by_character = strsplit(fasta,"")[[1]]
  covmat = matrix(0,length(by_character))
  Sequences = indind$Sequence
  Modifications = indind$Modifications
  Modifications=sub("N-Term",1,Modifications)
  match_ind=(str_locate_all(pattern = Sequences, fasta))
  modtype=regmatches(Modifications, gregexpr("(?<=\\().*?(?=\\))", Modifications, perl=T))
  
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
  stat_nr = as.matrix(table(totdf))
  
  totmat = matrix(0,nrow(covmat),nrow(stat_nr)+1)
  totmat[,1]=covmat

  for(i in 1:nrow(stat_nr)){
    totmat[as.numeric(colnames(stat_nr)),i+1]=stat_nr[i,]
  }

  totmat[,1]=totmat[,1]-rowSums(totmat[,2:ncol(totmat)])

  cols =c("darkblue","maroon","pink","brown","yellow","orange","red","blue","green","lightgrey")
  need= ncol(totmat)

  bp = barplot(t(totmat[,ncol(totmat):1]),col=cols[(length(cols)-need+1):length(cols)],border = NA)
  chars = rep("",nrow(totmat))
  isnot = which(rowSums(totmat[,2:ncol(totmat)])>0)
  chars[isnot]=by_character[isnot]
  mtext(side=1,at = bp,chars,cex=0.8)
  num=1:(length(by_character))
  #num = num-1
  num[!num%in%isnot]=""
  print(num)
  num = as.numeric(num)-1
  print(num)
  mtext(side=1,at = bp,num,cex=0.6,line=1,las=2)
  legend("topright",legend=c("no mod.",rownames(stat_nr)),fill=cols[length(cols):(length(cols)-need+1)])
  colnames(totmat)=c("no mod.",rownames(stat_nr))
  rownames(totmat)=paste0(by_character,0:(length(by_character)-1))
  return(totmat)
}

mm = ProteinPlot(fasta,indind) 

summaryFunction <- function(mm,indind){
  cov=sum(mm[,1]>0)/sum(mm[,1]>-1)
  counts=nrow(indind)
  modif=ncol(mm)-1
  modif.counts=sum(mm[,-1])
  mysum = data.frame(coverage=cov,nr.peptides=counts,modification_types=modif,modifications_detected=modif.counts)
  
}
summaryModif <- function(mm)
modonly = mm[rowSums(mm[,-1])>0,]
  
