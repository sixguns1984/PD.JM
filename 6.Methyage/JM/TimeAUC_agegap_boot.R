setwd("/public3/data/lijinxia/JM/JM/cg_batchcorrection/cg_Mvals/filter/weibull_AFT_aGH_scale_LEDD_PHS/intergration/methyage/QC_correct")
library(timeROC)
library(JM)
library(tidyverse)
library(parallel)
library(boot)
library(dplyr)

# data input
sampledata_PPMI <- read.csv("PPMI_jmdata_agegap_pred.csv", header = T)

# data prepared for timeROC
# PPMI
sampledata_PPMI <- sampledata_PPMI[order(sampledata_PPMI$Sample_Name,sampledata_PPMI$time),]

# bootstrap
timeROC.pre <- function(sampledata_PPMI, n){
  data_all <- sampledata_PPMI[,c(1:4,n)]
  data_all <- as.data.frame(data_all)
  data_BL <- subset(data_all, time == "0", select = c("Sample_Name", "status"))
  timeROC.boot <- function(data_BL, data_all, indices){
    data1 <- data_BL[indices,]
    data1 <- as.data.frame(data1)
    NDT <- data_all %>% as_tibble() %>% filter(.,Sample_Name %in% c(data1$Sample_Name))
    colnames(NDT)[5] <- "predict"
    ROC.pred <- timeROC(T = NDT$eventime, #piecewise-PH-aGH
                        delta = NDT$status,
                        marker = -NDT$predict,
                        cause = 1, weighting = "marginal",
                        times = seq(0,10,0.5),
                        iid = TRUE,
                        ROC = TRUE)
    return(ROC.pred$AUC)
  }
  set.seed(1234)
  ROC.boot <- boot(data = data_BL, statistic = timeROC.boot, R = 1000, strata = data_BL$status, data_all = data_all)
}
res <- mclapply(5:10, function(n)try(timeROC.pre(n, sampledata_PPMI = sampledata_PPMI), TRUE), mc.cores = 5)
save(res, file = "timeAUC_all_boot.RData")
