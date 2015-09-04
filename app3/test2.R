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

summaryFunction <- function(mm,indind){
  cov=sum(mm[,1]>0)/sum(mm[,1]>-1)
  counts=nrow(indind)
  modif=ncol(mm)-1
  modif.counts=sum(mm[,-1])
  mysum = data.frame(coverage=cov,nr.peptides=counts,modification_types=modif,modifications_detected=modif.counts)
  return(mysum)  
}

#summaryModif <- function(mm)
#modonly = mm[rowSums(mm[,-1])>0,]
  
