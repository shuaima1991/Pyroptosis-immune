rm(list = ls())
setwd("")
gene="riskScore"           
geneRT=read.table("gene.txt",sep="\t",header=F)          
files='CIBERSORT.result.txt'
data=read.table(files, header=T,sep="\t",check.names=F,row.names = 1)

##
CancerType=levels(as.factor(data$CancerType))
sameGenes=sameGenes=intersect(as.vector(geneRT[,1]),colnames(data))
outTab=matrix(1,32,22)

corTab=matrix(1,32,22)

for(i in 1:32){
  data1=data[data$CancerType==CancerType[i],]
  x=as.numeric(data1[,gene])
  
  for(j in 1:22){
    y=as.numeric(data1[,j])
    corT=cor.test(x,y)
    cor=corT$estimate
    pValue=corT$p.value
    
    outTab[i,j]=pValue
    corTab[i,j]=cor
  }}

colnames(outTab)=sameGenes
row.names(outTab)=CancerType
write.table(outTab,file="geneCor.pvalue.txt",sep="\t",row.names=T,quote=F)
colnames(corTab)=sameGenes
row.names(corTab)=CancerType
write.table(corTab,file="geneCor.cor.txt",sep="\t",row.names=T,quote=F)

