setwd("/public3/data/lijinxia/JM/JM/cg_batchcorrection/cg_Mvals/filter/weibull_AFT_aGH_scale_LEDD_PHS/intergration")
library(vroom)
sampledata <-vroom("/public3/data/lijinxia/JM/JM/cg_batchcorrection/cg_Mvals/filter/Sampledata015_LEDD_PHS_scale.csv")
sampledata <- as.data.frame(sampledata)
sampledata<-sampledata[order(sampledata$Sample.Event,decreasing = F),]

beta <- read.csv("/public3/data/lijinxia/JM/JM/cg_batchcorrection/cg_Mvals/filter/weibull_AFT_aGH_scale_LEDD_PHS/fitJOINT_bonf005.csv", header = T)
## bonferroni < 0.05
beta_bonf <- subset(beta, pvalue_bonferroni < 0.05)
beta_bonf05 <- beta_bonf[,2]
Mval_bonf <- sampledata[, match(beta_bonf$Name, colnames(sampledata))]
Mval_bonf <- as.matrix(Mval_bonf)

## FDR < 0.05
beta_fdr <- subset(beta, pvalue_FDR < 0.05)
beta_fdr05 <- beta_fdr[,2]
Mval_fdr <- sampledata[, match(beta_fdr$Name, colnames(sampledata))]
Mval_fdr <- as.matrix(Mval_fdr)

## beta from JM multiplies Mval at eatch point(bonferroni < 0.05)
Mval_bonf_inter <- as.data.frame(matrix(nrow = 552, ncol = 41))
names(Mval_bonf_inter) <- colnames(Mval_bonf)
for (n in 1:41){
  test <- Mval_bonf[,n] * beta_bonf05[n]
  Mval_bonf_inter[n] <- test
}
Mvalnew_bonf <- rowSums(Mval_bonf_inter)

## beta from JM multiplies Mval at eatch point(FDR < 0.05)
Mval_fdr_inter <- as.data.frame(matrix(nrow = 552, ncol = 34))
names(Mval_fdr_inter) <- colnames(Mval_fdr)
for (n in 1:34){
  test <- Mval_fdr[,n] * beta_fdr05[n]
  Mval_fdr_inter[n] <- test
}
Mvalnew_fdr <- rowSums(Mval_fdr_inter)


## merge Mval_inter and clinical message into one samplesheet
top2 <- Mval_bonf_inter[, match(beta_bonf$Name[1:2], colnames(Mval_bonf_inter))]

## build joint model for the new Mval
library(JM)
sampledata2$sex<-as.factor(sampledata2$sex)
sampledata2$time<-as.numeric(sampledata2$time)
sampledata2$LEDD<-scale(sampledata2$LEDD)
sampledata2<-sampledata2[order(sampledata2$Sample.Event,decreasing = F),]
sampledata2_BL<-subset(sampledata2,time==0)
colnames <- colnames(sampledata2)

# scaled the new Mval
sampledata2$Mvalnew_bonf <- scale(sampledata2$Mvalnew_bonf)
sampledata2$Mvalnew_fdr <- scale(sampledata2$Mvalnew_fdr)
i <- colnames[23]
fml <- as.formula( paste( i, "~", paste(c('age','time','sex','duration','CD8T','CD4T','NK','Bcell','Mono','Neu','LEDD'), collapse="+") ) )
fitLME <- lme(fml, data = sampledata2, random = ~1 | Sample_Name)
fitCOX <- coxph(Surv(eventime, status) ~ agediag + sex + educyrs + moca + GBA.APOE4.PHS, data = sampledata2_BL, x = TRUE)
fitJOINT <- jointModel(fitLME, fitCOX, timeVar = "time", method="weibull-AFT-aGH") # all of two < 0.05
res<-summary(fitJOINT)
res$'CoefTable-Event'

# predict the longitudinal part
newdata_PPMI <- sampledata2[,c('age','time','sex','duration','CD8T','CD4T','NK','Bcell','Mono','Neu','LEDD','agediag','educyrs','moca','GBA.APOE4.PHS','eventime','status','Sample_Name',i)]
newdata_PPMI <- newdata_PPMI[order(newdata_PPMI$Sample_Name,newdata_PPMI$time),]
lmpre_PPMI <- predict(fitJOINT, newdata_PPMI, type = c("Marginal"), interval = "confidence", return = TRUE)

# subset efficients of JM built from PPMI
coxcoef <- as.matrix(fixef(fitJOINT,process = "Event"))
lmpre_PPMI$sex <- as.numeric(lmpre_PPMI$sex)
coxdt2_PPMI <- as.matrix(lmpre_PPMI[,c('agediag','sex','educyrs', 'moca','GBA.APOE4.PHS','pred')])
jmpredict_PPMI <- as.matrix(coxdt2_PPMI) %*% as.matrix(coxcoef[2:7,]) + coxcoef[1,1]

sampledata3 <- cbind(sampledata[,1:22], Mvalnew_bonf, jmpredict_PPMI, top2)
names(sampledata3)[25:26] <- paste(colnames(top2), "jmpredict", sep = "_")
names(sampledata3)[24] <- "Mvalnew_bonf_jmpredict"
write.csv(sampledata3, "sampledata_Mvalnew.csv", row.names = F)
