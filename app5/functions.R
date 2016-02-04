#library("DESeq2")
load("genesAndcount.Rdata")
#genes= read.table("/Users/elin.axelsson/berger_group/user/elin.axelsson/analyses_to_change/cuffmerg_cuffquant_cuffnorm//Marchantia_csf2/analysis/genes.fpkm_table",row.names=1,header=T)
#samples = read.table("/Users/elin.axelsson/berger_group/user/elin.axelsson/analyses_to_change/cuffmerg_cuffquant_cuffnorm/Marchantia_csf2/analysis/samples.table",header=T,row.names=1)
#counts=read.table("/Users/elin.axelsson/berger_group/user/elin.axelsson/analyses_to_change/cuffmerg_cuffquant_cuffnorm/Marchantia_csf2/analysis/genes.count_table",row.names=1,header=T)
#s_names = sapply(strsplit(samples[,"file"],"/"),"[[",10)
#ind = (grep("r_",s_names))
#s_names[ind] = c("Mp_sperm1","Mp_sperm2","Mp_sperm3")
#colnames(genes)=s_names
#colnames(counts)=s_names
cols = data.frame(condition=c("10DAF","10DAF","15DAFembryo","15DAFembryo","20DAF1mmembryo","20DAF1mmembryo","20DAFSoftembryo","20DAFSoftembryo","female_3mm","female_3mm",
                              "Male_thalli","Male_thalli","Mp_female","Mp_thalli","Mp_sperm","Mp_sperm","Mp_sperm","Mp_sperm","Spore","Spore","Young_male_rep","Young_male_rep"))


notinc = list(
  c("Intercept","conditionYoung_male_rep"),
  c("Intercept","conditionMp_sperm"),
  c("Intercept","condition20DAF1mmembryo"),
  c("Intercept","condition15DAFembryo"),
  c("Intercept","condition20DAFSoftembryo"),
  c("Intercept","condition20DAF1mmembryo"),
  c("Intercept","conditionSpore"),
  c("Intercept","condition20DAFSoftembryo"),
  c("Intercept","condition15DAFembryo"),
  c("Intercept","condition10DAF")
)
ofInt = c("conditionMp_sperm","conditionYoung_male_rep","condition15DAFembryo","condition20DAF1mmembryo",
          "condition20DAF1mmembryo","condition20DAFSoftembryo","condition20DAFSoftembryo","conditionSpore","condition10DAF","condition15DAFembryo")



#####################################################

#counts = counts[1:500,]
#genes = genes[1:500,]

makeDataSet = function(data,outlier){
  torm = colnames(data)%in%outlier
  data = data[,!torm] 
  data
}
  
filterData = function(genes,counts,min=1){
  g_max = apply(genes,1,function(x) sum(x>min))
  keep = genes[which(g_max>1),]
  counts = counts[rownames(keep),]
  counts
}

runStats = function(counts,cols){
  dds <- DESeqDataSetFromMatrix(countData = round(counts),
                                colData = cols,
                                design = ~ condition)

  dds = estimateSizeFactors(dds)
  dds = estimateDispersions(dds)
  dds = nbinomLRT(dds,reduced=~1)
  res_anova =as.data.frame(results(dds))
  dds2= dds
  dds2 = nbinomWaldTest(dds2)    
  cpef = coefficients(dds2)
  return(list(cpef,res_anova))
}

generateStats = function(notinc,ofInt,cpef){
  test = cpef[,!colnames(cpef)%in%notinc]
  big = apply(test,1,function(x) max(x)==x[ofInt])
  mdist= apply(test,1,function(x) min(x[ofInt]-x[setdiff(colnames(test),ofInt)]))
  avdist= apply(test,1,function(x) mean(x[ofInt]-x[setdiff(colnames(test),ofInt)]))
  res = cbind(big,mdist,avdist)
  colnames(res) = paste(ofInt,colnames(res),sep="_")
  res
}

calcStatsAll = function(notinc,ofInt,cpef){
  statsList = list()
  for(i in 1:length(ofInt)){
    statsList[[i]]=generateStats(notinc[[i]],ofInt[i],cpef)
  }
  statsList
}

groupsIm = function(set,group,minFC,minavFC,st){
  if (group=="1"){
    i1 = 1
    i2 = 2
  }
  if (group=="2"){
    i1 = 3
    i2 = 4
  }
  if (group=="3"){
    i1 = 5
    i2 = 6
  }
  if (group=="4"){
    i1 = 7
    i2 = 8
  }
  if (group=="5"){
    i1 = 9
    i2 = 10
  }
  tg1 = names(which(st[[i1]][,1]>0 & st[[i1]][,2]>minFC & st[[i1]][,3]>minavFC))
  tg2 = names(which(st[[i2]][,1]>0 & st[[i2]][,2]>minFC & st[[i2]][,3]>minavFC))
  stat12 = length(intersect(tg1,tg2))
  stat1o2 = length(union(tg1,tg2))
  stat1 = length(setdiff(tg1,tg2))
  stat2 = length(setdiff(tg2,tg1))
  stats=c(stat1o2,stat12,stat1,stat2)
  
  if (set=="union"){
    tt = union(tg1,tg2)
  }  
  if (set=="intersect"){
    tt = intersect(tg1,tg2)
    
  }
  if (set=="setdiff1"){
    tt = setdiff(tg1,tg2)
    
  } 
  if (set=="setdiff2"){
    tt = setdiff(tg2,tg1)
    
  } 
  
  return(list(tt,stats))
}
  
groups = function(set,group,notinc,ofInt,cpef,minFC,minavFC){
  st=list()
  if (group=="1"){
    st[[1]] = generateStats(notinc[[1]],ofInt[1],cpef)
    st[[2]] = generateStats(notinc[[2]],ofInt[2],cpef)
  }
  if (group=="2"){
    st[[1]] = generateStats(notinc[[3]],ofInt[3],cpef)
    st[[2]] = generateStats(notinc[[4]],ofInt[4],cpef)
  }
  if (group=="3"){
    st[[1]] = generateStats(notinc[[5]],ofInt[5],cpef)
    st[[2]] = generateStats(notinc[[6]],ofInt[6],cpef)
  }
  if (group=="4"){
    st[[1]] = generateStats(notinc[[7]],ofInt[7],cpef)
    st[[2]] = generateStats(notinc[[8]],ofInt[8],cpef)
  }
  if (group=="5"){
    st[[1]] = generateStats(notinc[[9]],ofInt[9],cpef)
    st[[2]] = generateStats(notinc[[10]],ofInt[10],cpef)
  }
  
  tg1 = names(which(st[[1]][,1]>0 & st[[1]][,2]>minFC & st[[1]][,3]>minavFC))
  tg2 = names(which(st[[2]][,1]>0 & st[[2]][,2]>minFC & st[[2]][,3]>minavFC))
  stat12 = length(intersect(tg1,tg2))
  stat1o2 = length(union(tg1,tg2))
  stat1 = length(setdiff(tg1,tg2))
  stat2 = length(setdiff(tg2,tg1))
  stats=c(stat1o2,stat12,stat1,stat2)
  
  if (set=="union"){
    tt = union(tg1,tg2)
  }  
  if (set=="intersect"){
    tt = intersect(tg1,tg2)
   
  }
  if (set=="setdiff1"){
    tt = setdiff(tg1,tg2)
   
  } 
  if (set=="setdiff2"){
    tt = setdiff(tg2,tg1)
    
  } 
  
return(list(tt,stats))
}
  
####
load("PreCal.Rdata")

outls = sapply(prcal,"[[",3)
ds = lapply(prcal,"[[",1)
cl = lapply(prcal,"[[",2)
