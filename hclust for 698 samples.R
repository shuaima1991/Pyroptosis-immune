
library(sparcl)                                                                 
setwd("")                          
data=read.table("ssgseaOut.txt",sep="\t",header=T,check.names=F,row.names=1)    

#Delete normal samples
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[,group==0]

#hierarchical cluster
hc = hclust(dist(t(data)))

y=cutree(hc,2)               

write.table(y,file="cluster.txt",sep="\t",quote=F,col.names=F)

pdf(file="hclust.pdf",width=50,height=20)
ColorDendrogram(hc, y = y, labels = names(y), branchlength = 0.3,xlab=" ",sub=" ",main = " ")
dev.off()

