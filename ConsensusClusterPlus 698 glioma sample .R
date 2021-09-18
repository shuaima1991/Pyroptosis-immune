

setwd("D:/Rcode/����˼·/ϸ�������뿹��������/census ����")
#ȡ��R��ѧ�������Ĵ���
options(scipen = 1);
print(109000000);
rt=read.table("�����ļ�.txt",header=T,sep="\t",row.names = 1)
d=round(rt,2)
d=as.matrix(rt)

#���������оƬ������������һ���򵥵Ľ���normalization �����β�����λ�����Ļ�����Ȼ��ȡ�ڸ�����Ʒ����ܴ����Щgene����̽������������о������

mads=apply(d,1,mad)# mad(x) ������λ���� ���У�1��ȡd���ݵ���λ��

d=d[rev(order(mads))[1:13],]
#ȥ��ǰ5000������
d = sweep(d,1, apply(d,1,median,na.rm=T))
#���м�ȥ��λ����r������ʹ��sweep(x, MARGIN, STATS, FUN="-", ...) �Ծ���������㡣MARGINΪ1����ʾ�еķ����Ͻ������㣬
#Ϊ2��ʾ�еķ��������㡣STATS������Ĳ�����FUNΪ���㺯����Ĭ���Ǽ�����
colors <- colorRampPalette(c("white", "red"))(5)
library(ConsensusClusterPlus)
title=tempdir()
class(d)
results = ConsensusClusterPlus(d,maxK=6,reps=50,pItem=0.8,pFeature=1,
                               clusterAlg="hc",distance="pearson",plot="png",tmyPal =colors )

#������ĿK=2��3��4������������6�������س���������������80%������������β������ҵ��ȶ��ɿ���������ࡣ

#Ȼ��������Щ�����ǩ��������Ѱ�ҿ��Խ���������ı�ǩ���򡣿�������PAM������Ѱ��ǩ����

#results[[2]] is theresults result of k=2
results[[4]][["consensusMatrix"]][1:5,1:5]
results[[2]][["consensusTree"]]
#��k=����������������
results[[3]][["consensusClass"]]
c=write.table(results[[2]][["consensusClass"]],"3.txt",sep="\t")
icl = calcICL(results,title=title,plot="png")

icl[["clusterConsensus"]]
icl[["itemConsensus"]][1:5,]
#����result��list
c=results[[2]]$consensusClass
class(c)
c=as.data.frame(c)
b=write.table(icl,"2.txt",sep="\t")