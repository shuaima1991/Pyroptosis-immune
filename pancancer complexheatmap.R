setwd("D:/Rcode/文章思路/细胞焦亡与抗肿瘤免疫/泛癌与riskscore新/泛癌热图")
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
rt=read.table("9基因泛癌矩阵2.txt",header = T,row.names = 1,sep = "\t")
rt=as.matrix(rt)
rt=t(rt)

library(circlize)
rt=t(scale(t(rt)))
col_fun = colorRamp2(c(-1, 0, 1),c('blue', 'white', '#CD2626'))
#cluster_columns = FALSE使注释第一列按照聚类排列
Heatmap(rt, name = "mat", col = col_fun,show_column_names = FALSE,cluster_columns = FALSE,rect_gp = gpar(col = "white", lwd = 2))

