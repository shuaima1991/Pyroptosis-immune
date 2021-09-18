library(ggpubr)
pFilter=0.05
setwd("")                          
rt=read.table("inputdata.txt",sep="\t",header=T,check.names=F)    #fread data
rt=rt[!duplicated(rt$ID), ]
rownames(rt)=rt[,1]
rt=rt[,-1]
data=rt
data=log2(data)
#data=rt[rt[,"P-value"]<0.05,]
#data=data[,1:(ncol(rt)-3)]

Type=read.table("type.txt",sep="\t",check.names=F,header=F)
Type=Type[!duplicated(Type$V1), ]
rownames(Type)=Type[,1]
Type=Type[,-1]
Type=Type[row.names(data),]
colnames(Type)=c("cluster","Subtype")

outTab=data.frame()
data=cbind(data,Type)
for(i in colnames(data[,1:(ncol(data)-2)])){
  rt1=data[,c(i,"Subtype")]
  colnames(rt1)=c("expression","Subtype")
  ksTest<-kruskal.test(expression ~ Subtype, data = rt1)
  pValue=ksTest$p.value
    outTab=rbind(outTab,cbind(rt1,gene=i))
    print(pValue)
}
write.table(outTab,file="data.txt",sep="\t",row.names=F,quote=F)

#ggboxplot
data=read.table("data.txt",sep="\t",header=T,check.names=F)      
data$Subtype=factor(data$Subtype, levels=c("Normal","Glioma"))
p=ggboxplot(data, x="gene", y="expression", color = "Subtype",
            ylab="Expression value",
            xlab="",
            palette = c("red","blue") )
p=p+rotate_x_text(45)
pdf(file="nomal tumor.pdf",width=7,height=4)                         
p+stat_compare_means(aes(group=Subtype),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
dev.off()

