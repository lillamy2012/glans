library("DESeq2")

###############
## set up
###############
genes= read.table("/Users/elin.axelsson/berger_group/user/elin.axelsson/analyses_to_change/cuffmerg_cuffquant_cuffnorm//Marchantia_csf2/analysis/genes.fpkm_table",row.names=1,header=T)
samples = read.table("/Users/elin.axelsson/berger_group/user/elin.axelsson/analyses_to_change/cuffmerg_cuffquant_cuffnorm/Marchantia_csf2/analysis/samples.table",header=T,row.names=1)
counts=read.table("/Users/elin.axelsson/berger_group/user/elin.axelsson/analyses_to_change/cuffmerg_cuffquant_cuffnorm/Marchantia_csf2/analysis/genes.count_table",row.names=1,header=T)
s_names = sapply(strsplit(samples[,"file"],"/"),"[[",10)
ind = (grep("r_",s_names))
s_names[ind] = c("Mp_sperm1","Mp_sperm2","Mp_sperm3")
colnames(genes)=s_names
colnames(counts)=s_names
cols = data.frame(condition=c("10DAF","10DAF","15DAFembryo","15DAFembryo","20DAF1mmembryo","20DAF1mmembryo","20DAFSoftembryo","20DAFSoftembryo","female_3mm","female_3mm",
                              "Male_thalli","Male_thalli","Mp_female","Mp_thalli","Mp_sperm","Mp_sperm","Mp_sperm","Mp_sperm","Spore","Spore","Young_male_rep","Young_male_rep"))


o_counts=counts
o_genes = genes
o_cols = cols

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



#################
## functions
#################
makeDataSet = function(data,outlier){
  torm = colnames(data)%in%outlier
  data = data[,!torm] 
  data
}

calcStatsAll = function(notinc,ofInt,cpef){
  statsList = list()
  for(i in 1:length(ofInt)){
    statsList[[i]]=generateStats(notinc[[i]],ofInt[i],cpef)
  }
  statsList
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



##################################
### START 
##################################
sets = list(NULL,c("Mp_female","Mp_thalli","Mp_WT_sperm","Mp_sperm1"))
prcal = list()
i=0
for( outl in sets){
  i=i+1
  print(i)
  if (length(outl)>0){
    print("exclude")
    torm = colnames(counts)%in%outl   
    use_cols = cols[!torm,,drop=FALSE]
    use_counts = makeDataSet(counts,outl)
    use_genes = makeDataSet(genes,outl)
    use_counts = filterData(use_genes,use_counts,1)
  } else {
    print("default")
    use_counts = filterData(genes,counts,1)
    use_cols = cols
    use_genes = genes
  }
  print("runStats")
  r1 = runStats(use_counts,use_cols)
  stIn = r1[[1]][which(r1[[2]][,"padj"]<0.1),]
  print("calALL")
  cal1 = calcStatsAll(notinc,ofInt,stIn)
  print("add to list")
  prcal[[i]] = list(r1,cal1,outl)
}
save(prcal,file="PreCal.Rdata")
