setwd("/public3/data/lijinxia/JM/JM/cg_batchcorrection/cg_Mvals/filter/weibull_AFT_aGH_scale_LEDD_PHS/intergration/LMM/cutoff")
library(JM)
library(timeROC)
library(dplyr)
library(parallel)

sampledata <- read.csv("/public3/data/lijinxia/JM/JM/cg_batchcorrection/cg_Mvals/filter/weibull_AFT_aGH_scale_LEDD_PHS/intergration/sampledata_Mvalnew.csv", header = T)
sampledata$sex<-as.factor(sampledata$sex)
sampledata$time<-as.numeric(sampledata$time)
sampledata$LEDD<-scale(sampledata$LEDD)
sampledata$Mvalnew_bonf <- scale(sampledata$Mvalnew_bonf)
sampledata<-sampledata[order(sampledata$Sample.Event,decreasing = F),]
sampledata_BL<-subset(sampledata,time==0)
colnames <- colnames(sampledata)
 
## build joint model in PPMI
fitCOX <- coxph(Surv(eventime, status) ~ agediag + sex + educyrs + moca + GBA.APOE4.PHS, data = sampledata_BL, x = TRUE)
timeroc <- function(n, colnames, sampledata, sampledata_BL, fitCOX){
i <- colnames[n]
fml <- as.formula( paste( i, "~", paste(c('age','time','sex','duration','CD8T','CD4T','NK','Bcell','Mono','Neu','LEDD'), collapse="+") ) )
fitLME <- lme(fml,random=~1|Sample_Name,data = sampledata)
fitJOINT <- jointModel(fitLME, fitCOX, timeVar = "time", method="weibull-AFT-aGH")

#Calculates predicted values for the longitudinal part of a joint model
newdata_PPMI <- sampledata[,c(1:22,n)]
lmpre_PPMI <- predict(fitJOINT, newdata_PPMI, type = c("Marginal"), interval = "confidence", return = TRUE)

# subset efficients of JM built from PPMI
# PPMI
coxcoef <- as.matrix(fixef(fitJOINT,process = "Event")) #event coefficients
coxdt2_PPMI <- as.matrix(lmpre_PPMI[,c('agediag','sex','educyrs','moca','GBA.APOE4.PHS','pred')])
coxdt2_PPMI <- apply(coxdt2_PPMI, 2, as.numeric)
jmpredict_PPMI <- as.matrix(coxdt2_PPMI) %*% as.matrix(coxcoef[2:7,]) + coxcoef[1,]
lmedt_PPMI <- lmpre_PPMI[,c("Sample_Name","time","status","eventime")]
jmdata_PPMI <- cbind(lmedt_PPMI,jmpredict_PPMI)

# compute timeroc
ROC.pred.marginal <- timeROC(T = jmdata_PPMI$eventime, #piecewise-PH-aGH
                           delta = jmdata_PPMI$status,
                           marker = -jmdata_PPMI$jmpredict_PPMI,
                           cause = 1,weighting = "marginal",
                           times = c(6),
                           iid = TRUE,
                           ROC = TRUE)


# compute cutoff of timeroc from PPMI
jmdf_plot_PPMI <- data.frame(tpr = as.numeric(ROC.pred.marginal$TP),
                      fpr = as.numeric(ROC.pred.marginal$FP),
                      year = rep(c("0","6"),each = nrow(ROC.pred.marginal$TP)),
                      model = rep(c("joint model")))
jmdf_plot6_PPMI <- jmdf_plot_PPMI[which(jmdf_plot_PPMI$year == 6),]
jmdf_plot6_PPMI$youden <- as.numeric(jmdf_plot6_PPMI$tpr) - as.numeric(jmdf_plot6_PPMI$fpr)
cutoff6 <- jmdata_PPMI$jmpredict_PPMI[which.max(jmdf_plot6_PPMI$youden)]
cutoff6_coord <- jmdf_plot6_PPMI[which.max(jmdf_plot6_PPMI$youden),c(1:2)]
cutoff <- cbind(cutoff6, cutoff6_coord)
}
cutoff6 = mclapply(23:25, function(n) try(timeroc(n, colnames = colnames, sampledata = sampledata, sampledata_BL = sampledata_BL, fitCOX = fitCOX), TRUE), mc.cores = 4)
cutoff6 <- as.data.frame(cutoff6)
colnames(cutoff6) <- colnames[23:25]
test <- t(cutoff6)
test <- as.data.frame(test)
names(test)[1] <- "cutoff"
write.csv(test, file = "fitJOINT_Mvalnew_cutoff6.csv", row.names = T)


