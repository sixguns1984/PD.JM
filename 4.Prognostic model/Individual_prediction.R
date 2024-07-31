
setwd("D:/PPN_lab/PPN/硕士阶段/JM/JM_significant/filter/weibull_AFT_aGH_LEDD_PHS/individual_prediction")
library(JM)
library(tidyverse)

sampledata <- read.csv("sampledata_Mvalnew.csv", header = T)
sampledata$sex <- as.factor(sampledata$sex)
sampledata$time <- as.numeric(sampledata$time)
sampledata$Mvalnew_bonf <- scale(sampledata$Mvalnew_bonf)
sampledata$LEDD <- scale(sampledata$LEDD)

sampledata <- sampledata[order(sampledata$Sample_Name, sampledata$time),]

sampledata_BL <- subset(sampledata,time == 0)
colnames <- colnames(sampledata)

n <- which(colnames == 'Mvalnew_bonf') 

i <- colnames[n]

fml <- as.formula( paste( i, "~", paste(c('age','time','sex','duration','CD8T','CD4T','NK','Bcell','Mono','Neu','LEDD'), collapse="+") ) )
fitLME <- lme(fml,random=~1|Sample_Name,data = sampledata)
fitCOX <- coxph(Surv(eventime, status) ~ agediag + sex + educyrs + moca + GBA.APOE4.PHS, data = sampledata_BL, x = TRUE)
fitJOINT <- jointModel(fitLME, fitCOX, timeVar = "time", method="weibull-AFT-aGH")

sampledata_Mvalnew <- sampledata[,c('Sample_Name','eventime','status','age','time','sex','duration','CD8T','CD4T','NK','Bcell','Mono','Neu','LEDD','agediag','educyrs','moca','GBA.APOE4.PHS','Mvalnew_bonf')]
dataP51186 <- sampledata_Mvalnew[sampledata_Mvalnew$Sample_Name=='51186',]

#dynamic individual prediction
library(animation)
len_id <- nrow(dataP51186)
saveGIF({
  for(i in c(1:len_id)){
    sfit <- survfitJM(fitJOINT, newdata = dataP51186[1:i, ], idVar = "Sample_Name") 
    plot(sfit, estimator="mean", include.y = TRUE, conf.int=0.95, fill.area=TRUE, col.area="lightblue", main="Patient 51186")
    
  }
},ani.width = 400, ani.height=400)

#subset the importance (from junfeng)
#3125,3134,3312,50088,51186,51440
ND<-ND[ND$time<=3,]#依据此控制时间点
library(foreach)
myfun<-function(fitJOINT,sampledata_Mvalnew,i){
  a<-c('3108','3134','51186')
  ND <- sampledata_Mvalnew[sampledata_Mvalnew$Sample_Name == a[i], ]
  ND<-ND[ND$time <= 3,]
  ss <- survfitJM(fitJOINT, newdata = ND, idVar = "Sample_Name", M = 200, survTimes = seq(-1,8,1), last.time = 2)
  names(ss$summaries)<-'a'
  capns_sp<-data.frame(ss$summaries$a)
  capns_sp$pat<-a[i]
  capns_sp$group<-i
  return(capns_sp)
}
capns_sp<-foreach(i=1:3,.combine='rbind') %do% myfun(fitJOINT,sampledata_Mvalnew,i)
capns_sp$group <- factor(capns_sp$group)

#'3108','3134','51186'
longped<-function(fitJOINT,sampledata_Mvalnew,i){
  a<-c('3108','3134','51186')
  ND <- sampledata_Mvalnew[sampledata_Mvalnew$Sample_Name == a[i], ]
  ND<-ND[ND$time <= 3,]
  ss <- survfitJM(fitJOINT, newdata = ND, idVar = "Sample_Name", M = 200, survTimes = seq(-1,8,1), last.time = 2)
  names(ss$y)<-'a'
  names(ss$obs.times)<-'a'
  fit.y<-ss$y$a
  fit.t<-ss$obs.times$a
  lpr<-data.frame(fit.t,fit.y)
  lpr$pat<-a[i]
  lpr$group<-i
  return(lpr)
}
capns_lp<-foreach(i=1:3,.combine='rbind') %do% longped(fitJOINT,sampledata_Mvalnew,i)
capns_lp$group <- factor(capns_lp$group)

#dynamic individual prediction
library(ggplot2)
library(scales)
ggplot() +
  ylab('Epigenetic signature') + xlab(NULL) +
  geom_point(data = capns_lp, aes(x = fit.t, y = fit.y, colour = group), shape = 16, size = 2) +
  geom_smooth(data = capns_lp, aes(x = fit.t, y = fit.y, colour = group), se = FALSE, method = 'lm', size = 0.8) +
  coord_fixed(ratio = 0.8) +
  scale_y_continuous(expand = c(0, 0), limits = c(-3, 7), breaks = seq(-3, 7, by = 1), 
                     sec.axis = sec_axis( ~ rescale(., c(0, 100)),name = "Patient surviving free of MCI(%)")) +
  scale_x_continuous(expand = c(0, 0), limits = c(-0.3, 6), breaks = seq(0, 6, by = 2)) +
  scale_color_manual(label = c("Patient1","Patient2", "Patient3"),
                     values = c("#0073c2b2","#F39B7FB2",'#A73030B2')) +
  geom_vline(xintercept = 3, color = 'black', linetype = 2) +
  theme_classic() + #"Patient3",,'#A73030B2'
  theme(axis.text = element_text(size = 10))

ggplot(data = capns_sp, aes(x = times, y = Mean, colour = group, group = group)) + 
  ylab('Epigenetic signature') + xlab(NULL) +
  geom_smooth(method = "loess", span = 1, se = FALSE, size = 0.8) +
  # geom_ribbon(aes(ymin = Lower, ymax = Upper, fill=group), show.legend = FALSE, alpha = 0.5) +
  coord_fixed(ratio = 8) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1), breaks = seq(0, 1, by = 0.2),
                     sec.axis = sec_axis( ~ rescale(., c(0, 100)),name = "Patient surviving free of MCI(%)")) +
  scale_x_continuous(expand = c(0, 0), limits = c(-0.3, 6), breaks = seq(0, 6, by = 2)) +
  scale_color_manual(label = c("Patient1","Patient2", "Patient3")
                     ,values = c("#0073c2b2","#F39B7FB2",'#A73030B2')) +
  scale_fill_manual(values = c("#0073c2b2","#F39B7FB2",'#A73030B2')) + 
  theme_bw() + 
  geom_vline(xintercept = 3,color = 'black',linetype = 2) +
  theme_classic() #"Patient2","#F39B7FB2",
theme(axis.text = element_text(size = 10))


#roc plot at different time
ND <- sampledata[sampledata$Sample_Name == "3012", ]
roc <- rocJM(fitJOINT, dt = c(2, 4, 8), ND, idVar = "Sample_Name")
roc %>% plot()

