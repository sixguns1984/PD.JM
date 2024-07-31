
setwd("D:/PPN_lab/PPN/硕士阶段/JM/JM_significant/filter/weibull_AFT_aGH_LEDD_PHS/timeAUC")

load("timeAUC_bootstrap.RData")
load("timeAUC_PCClock_boot.RData")
#load("timeAUC_agegap_boot.RData")

## Wilcoxon rank sum test
pvalue_PPMI <- vector()
for (x in 4:19){
  wilcox_test <- wilcox.test(Troc_top$t[,x],Troc_inter$t[,x])
  pvalue_PPMI[x] <- wilcox_test$p.value
}
pvalue_PPMI

pvalue_Mvalnew <- vector()
for (x in 4:19){
  wilcox_test <- wilcox.test(Troc_inter$t[,x],Troc_inter_cox$t[,x])
  pvalue_Mvalnew[x] <- wilcox_test$p.value
}
pvalue_Mvalnew

# epigenetic signature and agegap
pvalue_Horv1 <- vector()
for (x in 4:19){
  wilcox_test <- wilcox.test(res[[1]]$t[,x],Troc_inter$t[,x])
  pvalue_Horv1[x] <- wilcox_test$p.value
}
pvalue_Horv1

pvalue_Horv2 <- vector()
for (x in 4:19){
  wilcox_test <- wilcox.test(res[[2]]$t[,x],Troc_inter$t[,x])
  pvalue_Horv2[x] <- wilcox_test$p.value
}
pvalue_Horv2

pvalue_Han <- vector()
for (x in 4:19){
  wilcox_test <- wilcox.test(res[[3]]$t[,x],Troc_inter$t[,x])
  pvalue_Han[x] <- wilcox_test$p.value
}
pvalue_Han

pvalue_Pheno <- vector()
for (x in 4:19){
  wilcox_test <- wilcox.test(res[[4]]$t[,x],Troc_inter$t[,x])
  pvalue_Pheno[x] <- wilcox_test$p.value
}
pvalue_Pheno

pvalue_Grim <- vector()
for (x in 4:19){
  wilcox_test <- wilcox.test(res[[4]]$t[,x],Troc_inter$t[,x])
  pvalue_Grim[x] <- wilcox_test$p.value
}
pvalue_Grim

# data prepared for timeauc line chart
PPMI_inter <- as.data.frame(Troc_inter$t0[4:18])
PPMI_inter$eventime <- substring(rownames(PPMI_inter), 3)
colnames(PPMI_inter) <- c("timeroc", "eventime")
PPMI_inter$group <- rep("Epigenetic signature", 15)

PPMI_Horv1 <- as.data.frame(res[[1]]$t0[4:18])
PPMI_Horv1$eventime <- substring(rownames(PPMI_Horv1), 3)
colnames(PPMI_Horv1) <- c("timeroc", "eventime")
PPMI_Horv1$group <- rep("Acceleration of PCHorvath1 age", 15)

PPMI_Horv2 <- as.data.frame(res[[2]]$t0[4:18])
PPMI_Horv2$eventime <- substring(rownames(PPMI_Horv2), 3)
colnames(PPMI_Horv2) <- c("timeroc", "eventime")
PPMI_Horv2$group <- rep("Acceleration of PCHorvath2 age", 15)

PPMI_Han2 <- as.data.frame(res[[3]]$t0[4:18])
PPMI_Han2$eventime <- substring(rownames(PPMI_Han2), 3)
colnames(PPMI_Han2) <- c("timeroc", "eventime")
PPMI_Han2$group <- rep("Acceleration of PCHannum age", 15)

PPMI_Pheno <- as.data.frame(res[[4]]$t0[4:18])
PPMI_Pheno$eventime <- substring(rownames(PPMI_Pheno), 3)
colnames(PPMI_Pheno) <- c("timeroc", "eventime")
PPMI_Pheno$group <- rep("Acceleration of PCPheno age", 15)

PPMI_Grim <- as.data.frame(res[[5]]$t0[4:18])
PPMI_Grim$eventime <- substring(rownames(PPMI_Grim), 3)
colnames(PPMI_Grim) <- c("timeroc", "eventime")
PPMI_Grim$group <- rep("Acceleration of PCGrim age", 15)

PPMI_auc <- rbind(PPMI_inter, PPMI_Horv1, PPMI_Horv2, PPMI_Han2, PPMI_Pheno, PPMI_Grim)

#time-dependent auc line chart
library(ggplot2)
library(timeROC)
PPMI_auc$group <- as.factor(PPMI_auc$group)
PPMI_auc$eventime <- as.numeric(PPMI_auc$eventime)

color <- c("#C59D94","#AFC7E8","#F09148","#427AB2","#DBDB8D","#FF9896")
#"#4682B4","#51B1B7","#E1C855","#E07B54"

ggplot(data = PPMI_auc, aes(x = eventime, y = timeroc, group = group, color = group)) +
  geom_point(size = 1, pch = 21, col = 'black') +
  coord_fixed(ratio = 15) +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) +
  geom_line(size = 1) +
  labs(x="Time since onset (years)", y="Time-dependent AUC") +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = c(0.27,0.3)) +
  scale_y_continuous(limits = c(0.5,1), breaks = seq(0.5,1,0.1), expand = c(0,0)) +
  scale_x_continuous(limits = c(1,8), breaks = seq(0,8,2), expand = c(0,0.2))


#第六年time-dependent ROC平滑曲线
library(ggplot2)
library(boot)
dync6 <- function(data, indices, colname){
  d <- data[indices,]
  ROC.pred.marginal<-timeROC(T=d$eventime, #piecewise-PH-aGH
                             delta=d$status,
                             marker=d$jmpredict,
                             cause=1,weighting="marginal",
                             times=c(6),
                             iid = TRUE,
                             ROC = TRUE)
  return(ROC.pred.marginal$AUC)
}
set.seed(1234)
Troc6 <- boot(data = jmdata, statistic = dync6, R = 5, stype = 'i', strata = jmdata$status)
result <- boot.ci(Troc6, type=c("perc", "bca"), index = 2)

coxdync6 <- function(data, indices){
  d <- data[indices,]
  ROC.coxpred.marginal<-timeROC(T=coxdata$eventime, 
                                delta=coxdata$status,
                                marker=coxdata$coxpredict,
                                cause=1,weighting="marginal",
                                times=c(6),
                                iid = TRUE,
                                ROC = TRUE)
  return(ROC.coxpred.marginal$AUC)
}
set.seed(1234)
cox_Troc6 <- boot(data = coxdata, statistic = coxdync6, R = 1000, stype = 'i', strata = coxdata$status)
coxresult <- boot.ci(cox_Troc6, type=c("perc", "bca"), index = 2)

wilcox.test()

# the 6th year
Troc6 <- Troc_inter$t[,13]
cox_Troc6 <- Troc_inter_cox$t[,13]
wilcox.test(Troc6, cox_Troc6) # p value < 2.2e-16
cutoff6 <- 0.792912180055735
cutoff6_coord <- c("0.9359","0")


library(ggplot2)
geom_point(aes(x = cutoff6_coord[2], y = cutoff6_coord[1]), size = 1, color = "black") 
geom_smooth(se=FALSE, size=1.2)  # 这就是平滑曲线的关键
theme_minimal(base_size = 14, base_family = "sans") 
annotate(geom = "point", x = cutoff6_coord$fpr, y = cutoff6_coord$tpr, colour = "red", size = 1.5) +
  annotate("text", x = cutoff6_coord$fpr, y = cutoff6_coord$tpr, label = paste("Cutoff: ", round(cutoff6,3)), color="black", size = 4 )
  
ggplot() +
  geom_line(data = df_plot, aes(x = fpr, y = tpr, color = model), size = 1.2) +
  geom_abline(slope = 1, intercept = 0, color = "grey10",linetype = 2, size = 1) +
  theme_classic() +
  coord_fixed(ratio = 1) +
  scale_color_manual(values = c("#E0367A","#0072B2"),
                     name = NULL, 
                     labels = c(paste0("AUC of joint model = ", round(AUC$jmAUC[14],2)), 
                                paste0("AUC of survival model = ",round(AUC$coxAUC[14],2)))) + 
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2), expand = c(0,0)) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2), expand = c(0,0)) +
  labs(x = "False positive rate (1 - specificity)", y = "True positive rate (sensitivity)") +
  theme(legend.position = c(0.7,0.15), 
        panel.border = element_blank(), #移除边界线
        panel.grid.major = element_blank(),  #移除网格线
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black")) 


#分组绘制jm的KM曲线
library(survival)
library(survminer)
jmdata$group <- ifelse(jmdata$jmpredict <= cutoff6, "low", "high")
jmdata$group <- as.factor(jmdata$group)
KM_fit <- survfit(Surv(eventime, status) ~ group, data = jmdata)
ggsurvplot(KM_fit, 
           data = jmdata, 
           pval = TRUE,
           pval.method = TRUE,
           palette = c("#E0367A","#0072B2"),
           risk.table = TRUE, 
           conf.int = TRUE)



