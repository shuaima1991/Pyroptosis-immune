setwd("")
library(maftools)
luad <- read.maf(maf="TCGA.GBM.LGG.maf")

# Extracting  correspondence"Tumor_Sample_Barcode" from clinical data
clin <- read.table("TCGA-GBM.LGG_phenotype .tsv", header=T, sep="\t")
clin.IS1 <- subset(clin, Type=="Lowimmunity  latepyroptosis")$Tumor_Sample_Barcode
clin.IS2 <- subset(clin, Type=="Lowimmunity  eralypyroptosis")$Tumor_Sample_Barcode
clin.IS3 <- subset(clin, Type=="Highimmunity  latepyroptosis")$Tumor_Sample_Barcode
clin.IS4 <- subset(clin, Type=="Highimmunity  eralypyroptosis")$Tumor_Sample_Barcode
# subsetMaf
luad.IS1 <- subsetMaf(maf=luad, tsb=clin.IS1, isTCGA=TRUE)
luad.IS2 <- subsetMaf(maf=luad, tsb=clin.IS2, isTCGA=TRUE)
luad.IS3 <- subsetMaf(maf=luad, tsb=clin.IS3, isTCGA=TRUE)
luad.IS4 <- subsetMaf(maf=luad, tsb=clin.IS4, isTCGA=TRUE)

#drive genes
shnsc.sigIS1 = oncodrive(maf = luad.IS1, AACol = 'HGVSp_Short', minMut = 5, pvalMethod = 'zscore')
plotOncodrive(res = shnsc.sigIS1, fdrCutOff = 0.01, useFraction = TRUE,labelSize = 1)
shnsc.sigIS2 = oncodrive(maf = luad.IS2, AACol = 'HGVSp_Short', minMut = 5, pvalMethod = 'zscore')
plotOncodrive(res = shnsc.sigIS2, fdrCutOff = 0.01, useFraction = TRUE,labelSize = 1)
shnsc.sigIS3 = oncodrive(maf = luad.IS3, AACol = 'HGVSp_Short', minMut = 5, pvalMethod = 'zscore')
plotOncodrive(res = shnsc.sigIS3, fdrCutOff = 0.01, useFraction = TRUE,labelSize = 1)
shnsc.sigIS4 = oncodrive(maf = luad.IS4, AACol = 'HGVSp_Short', minMut = 5, pvalMethod = 'zscore')
plotOncodrive(res = shnsc.sigIS4, fdrCutOff = 0.01, useFraction = TRUE,labelSize = 1)


#colour
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

oncoplot(luad.IS1, colors = vc_cols, top = 20)
oncoplot(luad.IS2, colors = vc_cols, top = 20)
oncoplot(luad.IS3, colors = vc_cols, top = 20)
oncoplot(luad.IS4, colors = vc_cols, top = 20)
output <- somaticInteractions(maf=luad.IS1, top=25, pvalue=c(0.05, 0.01))
output <- somaticInteractions(maf=luad.IS2, top=25, pvalue=c(0.05, 0.01))
output <- somaticInteractions(maf=luad.IS3, top=25, pvalue=c(0.05, 0.01))
output <- somaticInteractions(maf=luad.IS4, top=25, pvalue=c(0.05, 0.01))
output <- somaticInteractions(maf=luad,genes = c("CASP1","CASP3","CASP4","CASP5","CASP6","CASP7","CASP8","GSDMB","GSDMC","GSDMD","GSDME","GZMA","GZMB"), pvalue=c(0.05, 0.01))
#VAF
plotVaf(maf = luad.IS1,width = 100,height = 10)
plotVaf(maf = luad.IS2,width = 100,height = 10)
plotVaf(maf = luad.IS3,width = 100,height = 10)
plotVaf(maf = luad.IS4,width = 100,height = 10)
luad.vaf <- vafCompare(m1 = luad.IS2,m2 = luad.IS3)


# use mafCompare Compare differentially mutated genes
fvsm <- mafCompare(m1=luad.IS1, m2=luad.IS2, m3=luad.IS3, m4=luad.IS4, m1Name="IS1", m2Name="IS2",m3Name="IS3",m4Name="IS4", minMut=5)
# 结果保存到文件"female_vs_male.tsv"
write.table(fvsm$results, file="female_vs_male.tsv", quote=FALSE, row.names=FALSE, sep="\t")


#forest
fvsm1 <- mafCompare(m1=luad.IS1, m2=luad.IS2, m1Name="IS1", m2Name="IS2", minMut=5)
forestPlot(mafCompareRes=fvsm1, pVal=0.001, color=c("maroon", "royalblue"), geneFontSize=0.8)
fvsm2 <- mafCompare(m1=luad.IS3, m2=luad.IS4, m1Name="IS3", m2Name="IS4", minMut=5)
forestPlot(mafCompareRes=fvsm2, pVal=0.001, color=c("maroon", "royalblue"), geneFontSize=0.8)
fvsm3 <- mafCompare(m1=luad.IS1, m2=luad.IS4, m1Name="IS1", m2Name="IS4", minMut=5)
forestPlot(mafCompareRes=fvsm3, pVal=0.001, color=c("maroon", "royalblue"), geneFontSize=0.8)
fvsm4 <- mafCompare(m1=luad.IS2, m2=luad.IS3, m1Name="IS2", m2Name="IS3", minMut=5)
forestPlot(mafCompareRes=fvsm4, pVal=0.001, color=c("maroon", "royalblue"), geneFontSize=0.8)
fvsm5 <- mafCompare(m1=luad.IS1, m2=luad.IS3, m1Name="IS2", m2Name="IS3", minMut=5)
forestPlot(mafCompareRes=fvsm5, pVal=0.001, color=c("maroon", "royalblue"), geneFontSize=0.8)
#survival
mafSurvival(maf=luad.IS2, clinicalData = clin,genes=c("IDH1"), time="days_to_death.demographic", Status="vital_status.demographic", isTCGA=TRUE)
mafSurvival(maf=luad.IS3, clinicalData = clin,genes=c("PTEN"), time="days_to_death.demographic", Status="vital_status.demographic", isTCGA=TRUE)


lollipopPlot(maf = luad.IS2,maf = luad.IS3, gene = 'IDH1', AACol = 'HGVSp_Short', showMutationRate = TRUE)



lollipopPlot2(m1 = luad.IS2, m2 = luad.IS3, gene = "IDH1", AACol1 = "HGVSp_Short", AACol2 = "HGVSp_Short", m1_name = "Highimmunity  eralypyroptosis", m2_name = "Lowimmunity  latepyroptosis")

#output data
write.table(luad@variants.per.sample, file="TMB.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(luad@variant.type.summary, file="variant.type.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(luad@variant.classification.summary, file="variant.classification.txt", quote=FALSE, row.names=FALSE, sep="\t")

write.table(luad.IS2@variants.per.sample, file="IS2 TMB.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(luad.IS2@variant.type.summary, file="IS2 variant.type.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(luad.IS2@variant.classification.summary, file="IS2 variant.classification.txt", quote=FALSE, row.names=FALSE, sep="\t")

write.table(luad.IS3@variants.per.sample, file="IS3 TMB.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(luad.IS3@variant.type.summary, file="IS3 variant.type.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(luad.IS3@variant.classification.summary, file="IS3 variant.classification.txt", quote=FALSE, row.names=FALSE, sep="\t")

write.table(luad.IS4@variants.per.sample, file="IS4 TMB.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(luad.IS4@variant.type.summary, file="IS4 variant.type.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(luad.IS4@variant.classification.summary, file="IS4 variant.classification.txt", quote=FALSE, row.names=FALSE, sep="\t")

