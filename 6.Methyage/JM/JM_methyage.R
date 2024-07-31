
#!~/.conda/envs/R4.2.0/bin/R

setwd("/public3/data/lijinxia/JM/JM/cg_batchcorrection/cg_Mvals/filter/weibull_AFT_aGH_scale_LEDD_PHS/intergration/methyage/QC_correct")
library(JM)
sampledata<-read.csv("PPMI_sampledata_Mvalnew_methyage_correct_group.csv",header=T)

sampledata$sex<-as.factor(sampledata$sex)
sampledata$time<-as.numeric(sampledata$time)
sampledata$LEDD<-scale(sampledata$LEDD)
sampledata <- sampledata[order(sampledata$Sample_Name,sampledata$time),]
sampledata_BL<-subset(sampledata,time==0)

colnames<-colnames(sampledata)
colnames_age<-colnames[32:34]
fitCOX <- coxph(Surv(eventime, status) ~ agediag + sex + educyrs + moca + GBA.APOE4.PHS, data = sampledata_BL, x = TRUE)

pvalue<-vector('double',1)
abeta<-vector('double',1)
tbeta<-vector('double',1)
tpvalue<-vector('double',1)
Name<-c()
m=1

for(n in 1:3)try(
{
	i<-colnames_age[n]
	fml <- as.formula( paste( i, "~", paste(c('age','time','sex','duration','CD8T','CD4T','NK','Bcell','Mono','Neu',"LEDD"), collapse="+") ) )
	fitLME <- lme(fml,random=~1|Sample_Name,data = sampledata)
	# joint model fit
	trytemp<-try(fitJOINT <- jointModel(fitLME, fitCOX, timeVar = "time", method="weibull-AFT-aGH"),silent = TRUE)
	if ('try-error'%in% class(trytemp))
	  {
	pvalue[[m]]<-NA
 	abeta[[m]]<-NA
    	tpvalue[[m]]<-NA
   	tbeta[[m]]<-NA
    	Name[m]<-i
	m = m+1
    	next
  	} else {
	res<-summary(fitJOINT)
      	pvalue[[m]]<-res$'CoefTable-Event'['Assoct','p-value']
      	abeta[[m]]<-res$'CoefTable-Event'['Assoct','Value']
      	tbeta[[m]]<-res$'CoefTable-Long'['time','Value']
      	tpvalue[[m]]<-res$'CoefTable-Long'['time','p-value'] 
        Name[m]<-i 
	m = m+1
	}
})
fitJOINT<-data.frame(Name,abeta,pvalue,tbeta,tpvalue)
write.csv(fitJOINT,'weibull_AFT_aGH_scale_LEDD_PHS_agegap.csv',row.names=F)

