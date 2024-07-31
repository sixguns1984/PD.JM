setwd("/public3/data/lijinxia/JM/JM/cg_batchcorrection/cg_Mvals/filter/weibull_AFT_aGH_scale_LEDD_PHS/intergration/methyage/QC_correct")

library(timeROC)

# data input
sampledata<-read.csv("PPMI_sampledata_Mvalnew_methyage_correct_group.csv",header=T)
sampledata$sex<-as.factor(sampledata$sex)
sampledata$time<-as.numeric(sampledata$time)
sampledata$LEDD<-scale(sampledata$LEDD)
sampledata[,39:43] <- scale(sampledata[,39:43])
sampledata<-sampledata[order(sampledata$Sample.Event,decreasing = F),]
sampledata_BL<-subset(sampledata,time==0)

# build joint model for the Mvalnew_bonf
library(JM)
colnames <- colnames(sampledata)
i <- colnames[39:43]
agegap_pred = as.data.frame(matrix(nrow=552,ncol=5))

for (n in 1:5){
m <- i[n]
fitCOX <- coxph(Surv(eventime, status) ~ agediag + sex + educyrs + moca + GBA.APOE4.PHS, data = sampledata_BL, x = TRUE)
fml <- as.formula( paste( m, "~", paste(c('age','time','sex','duration','CD8T','CD4T','NK','Bcell','Mono','Neu','LEDD'), collapse="+") ) )
fitLME <- lme(fml,random=~1|Sample_Name,data = sampledata)
fitJOINT <- jointModel(fitLME, fitCOX, timeVar = "time", method="weibull-AFT-aGH")

# Calculates predicted values for the longitudinal part of a joint model
newdata <- sampledata[,c('age','time','sex','duration','CD8T','CD4T','NK','Bcell','Mono','Neu','LEDD','agediag','educyrs','moca','GBA.APOE4.PHS','Sample_Name','status','eventime',m)]
lmpre_inter <- predict(fitJOINT, newdata, type = c("Marginal"), interval = "confidence", return = TRUE)

# data prepared for timeROC of integrated joint model
coxcoef_inter <- as.matrix(fixef(fitJOINT,process = "Event")) #event coefficients
coxdt_inter2 <- as.matrix(lmpre_inter[,c('agediag','sex','educyrs', 'moca','GBA.APOE4.PHS', 'pred')])
coxdt_inter2 <- apply(coxdt_inter2, 2, as.numeric)
jmpredict_inter <- as.matrix(coxdt_inter2) %*% as.matrix(coxcoef_inter[2:7,]) + coxcoef_inter[1,1]
#lmedt_inter <- lmpre_inter[,c("Sample_Name","time","status","eventime")]
#jmdata_inter <- cbind(lmedt_inter,jmpredict_inter)
agegap_pred[,n] <- jmpredict_inter
}
colnames(agegap_pred) <- paste(i,"scale",sep="_")

sampledata2 <- sampledata[,c("Sample_Name","time","status","eventime")]
PPMI_jmdata_agegap_pred <- cbind(sampledata2, agegap_pred)
write.csv(PPMI_jmdata_agegap_pred, "PPMI_jmdata_agegap_pred.csv", row.names = F)
