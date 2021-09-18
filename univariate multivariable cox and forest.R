#Univariate
pFilter=0.05                                                   

setwd("")
library(survival)                                                 
library(UpSetR)
rt=read.table("inputdata.txt",header=T,sep="\t",check.names=F,row.names=1)       

outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
 cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
 coxSummary = summary(cox)
 coxP=coxSummary$coefficients[,"Pr(>|z|)"]
 outTab=rbind(outTab,
              cbind(id=i,
                    z=coxSummary$coefficients[,"z"],
                    HR=coxSummary$conf.int[,"exp(coef)"],
                    HR.95L=coxSummary$conf.int[,"lower .95"],
                    HR.95H=coxSummary$conf.int[,"upper .95"],
                    pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
              )
}
#Outputs all univariate results
outTab = outTab[is.na(outTab$pvalue)==FALSE,]
outTab=outTab[order(as.numeric(as.vector(outTab$pvalue))),]
write.table(outTab,file="uniCoxResult.txt",sep="\t",row.names=F,quote=F)

sigTab=outTab[as.numeric(as.vector(outTab$pvalue))<pFilter,]
write.table(sigTab,file="uniCoxResult.Sig.txt",sep="\t",row.names=F,quote=F)

sigGenes=c("futime","fustat")
sigGenes=c(sigGenes,as.vector(sigTab[,1]))
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="uniSigExp.txt",sep="\t",row.names=F,quote=F)


#Draw the upset
gene=sapply(strsplit(sigGenes,"\\|"),"[",1)
asType=sapply(strsplit(sigGenes,"\\|"),"[",3)
upsetList=list(AA=unique(gene[asType=="AA"]),
               AD=unique(gene[asType=="AD"]),
               AP=unique(gene[asType=="AP"]),
               AT=unique(gene[asType=="AT"]),
               ES=unique(gene[asType=="ES"]),
               ME=unique(gene[asType=="ME"]),
               RI=unique(gene[asType=="RI"]) )
upsetData=fromList(upsetList)

pdf(file="uniCoxUpset.pdf",onefile = FALSE,width=8,height=5)              
upset(upsetData,
      nsets = 7,                                    
      order.by = "freq",                            
      show.numbers = "yes",                        
      number.angles = 20,                           
      point.size = 1.5,                            
      matrix.color="red",                           
      line.size = 0.8,                              
      mainbar.y.label = "Gene Intersections",
      sets.x.label = "Set Size")
dev.off()


#multivariable
setwd("") 

library(survival)
library(survminer)

rt=read.table("inputdata",header=T,sep="\t",check.names=F,row.names=1)
#rt[,"VCAN"]=log2(rt[,"VCAN"]+1)

rt$Survival_month=rt$Survival_month+0.04
multiCox=coxph(Surv(rt$Survival_month, rt$status) ~ ., data = rt)
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
   HR=multiCoxSum$conf.int[,"exp(coef)"],
   HR.95L=multiCoxSum$conf.int[,"lower .95"],
   HR.95H=multiCoxSum$conf.int[,"upper .95"],
   pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox2.xls",sep="\t",row.names=F,quote=F)

pdf(file="forest.pdf",
    width = 7,             
    height = 6,            
)
ggforest(multiCox,
         main = "Hazard ratio",
         cpositions = c(0.02,0.22, 0.4), 
         fontsize = 0.7, 
         refLabel = "reference", 
         noDigits = 2)
dev.off()

#Forest
setwd("") 
library(forestplot)

bioForest=function(coxFile=null,forestFile=null,forestCol=null){
   
   rt=read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
   data=as.matrix(rt)
   HR=data[,1:3]
   hr=sprintf("%.3f",HR[,"HR"])
   hrLow=sprintf("%.3f",HR[,"HR.95L"])
   hrHigh=sprintf("%.3f",HR[,"HR.95H"])
   pVal=data[,"pvalue"]
   pVal=ifelse(pVal<0.05, "<0.05", sprintf("%.3f", pVal))
   clrs <- fpColors(box=forestCol,line="darkblue", summary="royalblue")      #定义颜色
   tabletext <- 
      list(c(NA, rownames(HR)),
           append("pvalue", pVal),
           append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )   #定义图片文字
   pdf(file=forestFile,width = 9,height = 10,onefile = FALSE)
   forestplot(tabletext, 
              rbind(rep(NA, 3), HR),
              col=clrs,
              graphwidth=unit(50, "mm"),
              xlog=T,
              lwd.ci=4,
              boxsize=0.6,
              xlab="Hazard ratio",
              txt_gp=fpTxtGp(ticks=gpar(cex=1.1),xlab=gpar(cex = 1.25))
   )
   dev.off()
}
############Draw the forest map ############

bioForest(coxFile="multiCox2.txt",forestFile="forest3.pdf",forestCol="red")

