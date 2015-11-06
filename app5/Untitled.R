genes= read.table("../../../Marchantia_csf2/analysis/genes.fpkm_table",row.names=1,header=T)
samples = read.table("../../../Marchantia_csf2/analysis/samples.table",header=T,row.names=1)
counts=read.table("../../../Marchantia_csf2/analysis/genes.count_table",row.names=1,header=T)
## rpkm to TPM 
#tpm = apply(genes,2,function(x) y=10^6*x/(sum(x)))
s_names = sapply(strsplit(samples[,"file"],"/"),"[[",10)
#colnames(tpm)=s_names
colnames(genes)=s_names
g_max = apply(genes,1,function(x) sum(x>1))
keep = genes[which(g_max>1),]

cols = data.frame(condition=c("10DAF","10DAF","15DAFembryo","15DAFembryo","20DAF1mmembryo","20DAF1mmembryo","20DAFSoftembryo","20DAFSoftembryo","female_3mm","female_3mm",
                               "Male_thalli","Male_thalli","Mp_female","Mp_thalli","Mp_sperm","Mp_sperm","Mp_sperm","Mp_sperm","Spore","Spore","Young_male_rep","Young_male_rep"))
counts = counts[rownames(keep),]
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


generateStats = function(notinc,ofInt,cpef){
  test = cpef[,!colnames(cpef)%in%notinc]
  big = apply(test,1,function(x) max(x)==x[ofInt])
  mdist= apply(test,1,function(x) min(x[ofInt]-x[setdiff(colnames(test),ofInt)]))
  avdist= apply(test,1,function(x) mean(x[ofInt]-x[setdiff(colnames(test),ofInt)]))
  res = cbind(big,mdist,avdist)
  colnames(res) = paste(ofInt,colnames(res),sep="_")
  res
}

notinc = list(
  c("Intercept","conditionYoung_male_rep"),
  c("Intercept","conditionMp_sperm"),
  c("Intercept","condition20DAF1mmembryo"),
  c("Intercept","condition10DAF"),
  c("Intercept","condition20DAFSoftembryo"),
  c("Intercept","condition20DAF1mmembryo"),
  c("Intercept","conditionSpore"),
  c("Intercept","condition20DAFSoftembryo")
)
ofInt = c("conditionMp_sperm","conditionYoung_male_rep","condition10DAF","condition20DAF1mmembryo","condition20DAF1mmembryo","condition20DAFSoftembryo","condition20DAFSoftembryo","conditionSpore")

statsList = list()
for(i in 1:length(ofInt)){
  statsList[[i]]=generateStats(notinc[[i]],ofInt[i],cpef)
}

sttot = do.call("cbind",statsList)
tot = cbind(keep,res_anova$padj,sttot)
colnames(tot)[duplicated(colnames(tot))] = paste("v2",colnames(tot)[duplicated(colnames(tot))],sep="_")
tot$group= 0
gr1 = which(tot$conditionMp_sperm_big+tot$conditionYoung_male_rep_big>0)
gr2 = which(tot$condition10DAF_big+tot$condition20DAF1mmembryo_big>0)
gr3 = which(tot$condition20DAFSoftembryo_big+tot$v2_condition20DAF1mmembryo_big>0)
gr4 = which(tot$v2_condition20DAFSoftembryo_big+tot$conditionSpore_big>0)
tot$group[gr1]=1
tot$group[gr2]=2
tot$group[gr3]=3
tot$group[gr4]=4
tot = tot[,c("young_male_reprod_r1", "young_male_reprod_r2", "Mp_WT_sperm", "r_1","r_2" ,"r_3","10DAF_r1", "10DAF_r2", "15DAFembry1" ,"15DAFembry2",
             "20DAF1mmembry1","20DAF1mmembry2","20DAFSoftembry1","20DAFSoftembry2" ,"spore1","spore2","female_3mm_r1", "female_3mm_r2","Mp_female", 
             "Male_thalli_r1", "Male_thalli_r2", "Mp_thalli","res_anova$padj",
             "conditionMp_sperm_big","conditionMp_sperm_mdist","conditionMp_sperm_avdist", "conditionYoung_male_rep_big" ,"conditionYoung_male_rep_mdist",
          "conditionYoung_male_rep_avdist", "condition10DAF_big", "condition10DAF_mdist", "condition10DAF_avdist", "condition20DAF1mmembryo_big", "condition20DAF1mmembryo_mdist"
          ,"condition20DAF1mmembryo_avdist", "v2_condition20DAF1mmembryo_big", "v2_condition20DAF1mmembryo_mdist","v2_condition20DAF1mmembryo_avdist",
"v2_condition20DAFSoftembryo_mdist","v2_condition20DAFSoftembryo_avdist", "conditionSpore_big", "conditionSpore_mdist", "conditionSpore_avdist", "group")]
colnames(tot)[4:6] = c("Mp_sperm1","Mp_sperm2","Mp_sperm3")

save(tot,file="table.Rdata")

#toUse = test[which(to$big==1 & test$`res_anova$padj`<0.05 & test$mdist>2),]
