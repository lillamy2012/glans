load("data/rpkm_Mp.rdata")
conds =as.factor(c(rep("10DAF",2),rep("15DAFemb",2),rep("20DAF1mm",2),rep("20DAFSoft",2),rep("female_3mm",2),
rep("male_thalli",2),"female","thalli","sperm",rep("spore",2),rep("ymalepre",2)))

quantileData <- function(data,n=10){
    qq=quantile(data[which(data>0)],seq(0,1,length.out = n))
    cdata= as.numeric(cut(data,c(0,qq[-1]),include.lowest=F))
    cdata[is.na(cdata)]=0
    cdata
}
meanSdData <- function(data=a[,-1],set=conds){
    m=matrix(NA,nrow=nrow(a),ncol=nlevels(conds))
    colnames(m)=levels(conds)
    s = m
    for (i in levels(conds)){
        m[,i]= rowMeans(data[,which(conds==i),drop=F])
        s[,i] = apply(data[,which(conds==i),drop=F],1,sd)
    }
    tot=list(m,s)
}
qdata = apply(a[,-1],2,quantileData)
meansMatS = meanSdData(a[,-1],set=conds)
meansMat = data.frame(loc=a[,1],meansMatS[[1]])
sMat = data.frame(loc=a[,1],meansMatS[[2]])
save(a,qdata,meansMat,sMat,file="data/app1_data.rdata")
