
setwd("")
library(ggplot2)
library(ggpubr)
library(ggExtra)
exp=read.table("inputdata.txt", header=T,sep="\t",check.names=F)
x=exp$riskScore
y=exp$TSI
corT=cor.test(x,y,method="spearman")
cor=corT$estimate
pValue=corT$p.value
df1=as.data.frame(cbind(x,y))
df1=log2(df1)
p1=ggplot(df1, aes(x, y)) + 
  xlab(x)+ylab(y)+
  ggtitle(paste0("Cancer: ","x"))+theme(title=element_text(size=10))+
  geom_point()+ geom_smooth(method="lm") + theme_bw()+
  stat_cor(method = 'spearman', aes(x =x, y =y))
p2=ggMarginal(p1, type = "density", xparams = list(fill = "orange"),yparams = list(fill = "blue"))
print(p2)
