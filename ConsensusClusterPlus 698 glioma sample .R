setwd("")
options(scipen = 1);
print(109000000);
rt=read.table("输入文件.txt",header=T,sep="\t",row.names = 1)
d=round(rt,2)
d=as.matrix(rt)
mads=apply(d,1,mad)# mad(x) 

d=d[rev(order(mads))[1:13],]

d = sweep(d,1, apply(d,1,median,na.rm=T))

colors <- colorRampPalette(c("white", "red"))(5)
library(ConsensusClusterPlus)
title=tempdir()
class(d)
results = ConsensusClusterPlus(d,maxK=6,reps=50,pItem=0.8,pFeature=1,
                               clusterAlg="hc",distance="pearson",plot="png",tmyPal =colors )

#results[[2]] is theresults result of k=2
results[[4]][["consensusMatrix"]][1:5,1:5]
results[[2]][["consensusTree"]]
results[[3]][["consensusClass"]]
c=write.table(results[[2]][["consensusClass"]],"3.txt",sep="\t")
icl = calcICL(results,title=title,plot="png")

icl[["clusterConsensus"]]
icl[["itemConsensus"]][1:5,]
c=results[[2]]$consensusClass
class(c)
c=as.data.frame(c)
b=write.table(icl,"2.txt",sep="\t")
