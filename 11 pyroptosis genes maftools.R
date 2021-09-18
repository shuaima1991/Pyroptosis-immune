setwd("")
library(maftools)
luad <- read.maf(maf="TCGA.GBM.LGG合并�?.maf",clinicalData = "TCGA-GBM.LGG_phenotype .tsv")

#
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
print(vc_cols)
#>   Frame_Shift_Del Missense_Mutation Nonsense_Mutation         Multi_Hit 
#>         "#A6CEE3"         "#1F78B4"         "#B2DF8A"         "#33A02C" 
#>   Frame_Shift_Ins      In_Frame_Ins       Splice_Site      In_Frame_Del 
#>         "#FB9A99"         "#E31A1C"         "#FDBF6F"         "#FF7F00"
oncoplot(luad, colors = vc_cols, genes = c("CASP1","CASP3","CASP4","CASP5","CASP11","CASP8","GSDMB","GSMDC","GSDMD","GSDME","GZMA","GZAMB"))
