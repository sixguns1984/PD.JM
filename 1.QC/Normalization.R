#!/home/lijinxia/bin/R

library(minfi)

#set up a path to the data directory
dataDirectory = '/public/labdata/lijinxia/JM/ppmi_140/'

#read in the sample sheet for the experiment
targets <- read.metharray.sheet(dataDirectory, pattern="Sample_Sheet_PPMI_140_data.csv")

#read in the raw data from the IDAT files
RGset <- read.metharray.exp(targets=targets,force = TRUE)
RGset

#give the samples descriptive names
targets$Sample_ID <- paste(targets$Sample_Group,targets$Sample_Name,sep='-')
sampleNames(RGset) <- targets$Sample_ID

#calculate the detection p-values
detP <- detectionP(RGset)

#examine mean detection p-values across all samples to identify any failed samples
library(RColorBrewer)
pal <- brewer.pal(8,"Dark2") 
barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2,cex.names=0.8,ylim = c(0,0.01), ylab="Mean detection p-values") 
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal,bg="white", bty = "n", x.intersp = 0.6, text.width = 1.4, y.intersp = 0.8)

#remove samples with low quality£¨p>0.01£©
keep <- colMeans(detP) < 0.01 
RGset1 <- RGset[,keep]
targets1 <- targets[keep,]
detP1 <- detP[,keep]

Mset = preprocessRaw(RGset1)
Gset = mapToGenome(Mset)  
annotation = getAnnotation(Gset)

##predict sex
pre_sex <- getSex(Gset)
Sample_Sheet_PPMI_140_data = read.csv("/public/labdata/lijinxia/JM/ppmi_140/Sample_Sheet_PPMI_140_data.csv",header = TRUE,check.names = F)
Sample_Sheet_PPMI_140_data$presex <- pre_sex$predictedSex
Sample_Sheet_PPMI_140_data2<-subset(Sample_Sheet_PPMI_140_data,Sample_Sheet_PPMI_140_data$SEX==Sample_Sheet_PPMI_140_data$presex)

RGset2 <- RGset1[,c(Sample_Sheet_PPMI_140_data2$ID_Position)]
Mset1 <- Mset[,c(Sample_Sheet_PPMI_140_data3$ID_Position)]

pdf(file = "presex.pdf")
plotSex(getSex(Gset1, cutoff = -2))
dev.off()

##detect abnormal samples£¨uMeth/mMeth < 10.5)
qc <- getQC(Mset1)
pdf (file = "QC.pdf")
plotQC(qc, badSampleCutoff = 10.5)
dev.off()
Mset2 <- addQC(Mset1, qc = qc)

keep1 <- (colnames(RGset2) %in% Mset2$ID_Position)
RGset3 <- RGset2[,keep1]

##MDSplot
pd <- pData(RGset3)
names <- pd$Sample_Name
groups <- pd$Sample_Group
pdf(file = "MDSplotraw.pdf")
mdsPlot(RGset3,numPositions = 1000, sampNames=names, sampGroups=groups)
dev.off()


##probes with a detection P value >0.01 in 20% or more of samples were identified
detP1 = detectionP(RGset3)
keep2 = apply(detP1,1,function(x) {all(x < 0.01)})
RGset4 = RGset3[keep2,]


##normalization
Mset.noob <- preprocessNoob(RGset4)
Beta.noob = getBeta(Mset.noob1)
write.csv(Beta.noob,'/public/labdata/lijinxia/JM/QC/sixth_progression/methy_beta.csv')

#Test the effect of normalization
library(RColorBrewer)
pd1 <- pData(RGset4)
pdf(file = "Raw beta1.pdf")
densityPlot(RGset4, sampGroups=pd1$Sample_Group,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(pd1$Sample_Group)), x.intersp = 0.6, cex = 0.6, text.col=brewer.pal(8,"Dark2"))
dev.off()
pd3 <- pData(Mset.noob)
pdf(file = "Normalized beta.pdf")
densityPlot(Mset.noob, sampGroups=pd3$Sample_Group,main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(pd2$Sample_Group)), x.intersp = 0.6, cex = 0.6, text.col=brewer.pal(8,"Dark2"))
dev.off()


#filter data
#Filter probes with low confidence
Gset.noob = mapToGenome(Mset.noob1)
detP2 <- detP1[match(featureNames(Gset.noob),rownames(detP1)),] 
keep3 <- rowSums(detP2 < 0.01) == ncol(Gset.noob) 
table(keep3)
GSetSqFlt <- Gset.noob[keep3,]

#Filter probes on sex chromosome
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
keep4 <- !(featureNames(GSetSqFlt) %in% ann850k$Name[ann850k$chr %in% c("chrX","chrY")])
table(keep4)
GSetSqFlt1 <- GSetSqFlt[keep4,]

#Filter probes related with snp
GSetSqFlt2 <- dropLociWithSnps(GSetSqFlt1)
GSetSqFlt2

#Filter probes mapped to multiple locations
library(maxprobes)
xloci <- maxprobes::xreactive_probes(array_type = "EPIC")
length(xloci)
GSetSqFlt3 <- maxprobes::dropXreactiveLoci(GSetSqFlt2)
GSetSqFlt3


#MDS plot
pd3 <- pData(GSetSqFlt3)
names <- pd3$Sample_Name
groups <- pd3$Sample_Group
pdf(file = "MDSplotnormalized.pdf")
mdsPlot(getBeta(GSetSqFlt3),numPositions = 1000, sampNames=names, sampGroups=groups)
dev.off()

##Cell component
library(FlowSorted.Blood.EPIC)
countsEPIC <- estimateCellCounts2(RGset4, compositeCellType = "Blood", processMethod = "preprocessNoob", cellTypes = c("CD8T","CD4T","NK","Bcell","Mono", "Neu"),referencePlatform = c("IlluminaHumanMethylationEPIC"))

countsEPIC <- estimateCellCounts2(RGset4,localHub=TRUE, compositeCellType = "Blood", processMethod = "preprocessNoob", cellTypes = c("CD8T","CD4T","NK","Bcell","Mono", "Neu"),referencePlatform = c("IlluminaHumanMethylationEPIC"))


#Compute methylation level
mVals <- getM(GSetSqFlt3)
bVals <- getBeta(GSetSqFlt3)
write.csv(bVals, file = "beta values.csv")
write.csv(mVals, file = "M values.csv")

#Density plot
pdf(file = "beta values.pdf")
densityPlot(bVals, sampGroups=pd2$Sample_Group, main="Beta values", legend=FALSE, xlab="Beta values")
legend("top", legend = levels(factor(pd2$Sample_Group)), x.intersp = 0.6, cex = 0.6, text.col=brewer.pal(8,"Dark2"))
dev.off()
pdf(file = "M values.pdf")
densityPlot(mVals, sampGroups=pd2$Sample_Group, main="M-values", legend=FALSE, xlab="M values")
legend("top", legend = levels(factor(pd2$Sample_Group)), x.intersp = 0.6, cex = 0.6, text.col=brewer.pal(8,"Dark2"))
dev.off()
