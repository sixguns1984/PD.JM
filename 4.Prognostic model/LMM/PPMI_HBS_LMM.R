setwd("/public3/data/lijinxia/JM/JM/cg_batchcorrection/cg_Mvals/filter/weibull_AFT_aGH_scale_LEDD_PHS/intergration/LMM")
library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest) #calculate p value

## data of cutoff6 impute
cutoff6 <- read.csv("/public3/data/lijinxia/JM/JM/cg_batchcorrection/cg_Mvals/filter/weibull_AFT_aGH_scale_LEDD_PHS/intergration/LMM/cutoff/fitJOINT_Mvalnew_cutoff6.csv", header = T)
colnames(cutoff6) <- c("Name", "cutoff")

## data input from PPMI
sampledata_PPMI <- read.csv("/public3/data/lijinxia/JM/JM/cg_batchcorrection/cg_Mvals/filter/weibull_AFT_aGH_scale_LEDD_PHS/intergration/sampledata_Mvalnew.csv", header = T)

## data input from HBS
sampledata_HBS <- read.csv("sampledata_Mvalnew_HBS.csv", header = T)
sampledata_HBS <- rename(sampledata_HBS, Sample.Event = Sample_Name, Sample_Name = ID)

### divide corhort into different groups according to the cutoff6 at maximum follow-up time
# PPMI top2 CpG
group_PPMI <- sampledata_PPMI[order(sampledata_PPMI$Sample_Name, sampledata_PPMI$time),]
group_PPMI <- group_PPMI %>% group_by(Sample_Name) %>% slice(which.max(time))
group_PPMI$group1 <- ifelse(group_PPMI$cg08477640_jmpredict <= cutoff6[2,2], "low", "high")
group_PPMI$group2 <- ifelse(group_PPMI$cg13241411_jmpredict<= cutoff6[3,2], "low", "high")
group_PPMI <- group_PPMI[, c("Sample_Name", "group1", "group2")]
sampledata_PPMI2 <- merge(sampledata_PPMI, group_PPMI, by = "Sample_Name", all.x = T)

sampledata_PPMI2$age <- 10^(sampledata_PPMI2$age)
sampledata_PPMI2$sex <- as.factor(sampledata_PPMI2$sex)
sampledata_PPMI2$Sample_Name <- as.factor(sampledata_PPMI2$Sample_Name)
sampledata_PPMI2$group1 <- as.factor(sampledata_PPMI2$group1)
sampledata_PPMI2$group2 <- as.factor(sampledata_PPMI2$group2)
sampledata_PPMI2 <- sampledata_PPMI2[order(sampledata_PPMI2$Sample.Event, decreasing = F),]


# PPMI integrated Mval
group_PPMI2 = sampledata_PPMI[order(sampledata_PPMI$Sample_Name, sampledata_PPMI$time),]
group_PPMI2 <- group_PPMI2 %>% group_by(Sample_Name) %>% slice(which.max(time))
group_PPMI2$group <- ifelse(group_PPMI2$Mvalnew_bonf_jmpredict <= cutoff6[1,2], "low", "high")
group_PPMI2 <- group_PPMI2[, c("Sample_Name", "group")]
sampledata_PPMI3 <- merge(sampledata_PPMI, group_PPMI2, by = "Sample_Name", all.x = T)

sampledata_PPMI3$age <- 10^(sampledata_PPMI3$age)
sampledata_PPMI3$sex <- as.factor(sampledata_PPMI3$sex)
sampledata_PPMI3$Sample_Name <- as.factor(sampledata_PPMI3$Sample_Name)
sampledata_PPMI3$group <- as.factor(sampledata_PPMI3$group)
sampledata_PPMI3 <- sampledata_PPMI3[order(sampledata_PPMI3$Sample.Event, decreasing = F),]

# HBS top2 CpG
group_HBS <- sampledata_HBS[order(sampledata_HBS$Sample_Name, sampledata_HBS$time),]
group_HBS <- group_HBS %>% group_by(Sample_Name) %>% slice(which.max(time))
group_HBS$group1 <- ifelse(group_HBS$cg08477640_jmpredict <= cutoff6[2,2], "low", "high")
group_HBS$group2 <- ifelse(group_HBS$cg13241411_jmpredict <= cutoff6[3,2], "low", "high")
group_HBS <- group_HBS[, c("Sample_Name", "group1", "group2")]
sampledata_HBS2 <- merge(sampledata_HBS, group_HBS, by = "Sample_Name", all.x = T)

sampledata_HBS3 <- sampledata_HBS2[,c("Sample_Name","time","status","duration","MMSE","GBA.APOE4.PHS","sex","age","educyrs","group1","group2")]
sampledata_HBS3$sex <- ifelse(sampledata_HBS3$sex == "M", "1", "2")
sampledata_HBS3$sex <- as.factor(sampledata_HBS3$sex)
sampledata_HBS3$Sample_Name <- as.factor(sampledata_HBS3$Sample_Name)
sampledata_HBS3$group1 <- as.factor(sampledata_HBS3$group1)
sampledata_HBS3$group2 <- as.factor(sampledata_HBS3$group2)
sampledata_HBS3 <- sampledata_HBS3[order(sampledata_HBS3$Sample_Name, sampledata_HBS3$time),]

# HBS integrated Mval
group_HBS2 <- sampledata_HBS[order(sampledata_HBS$Sample_Name, sampledata_HBS$time),]
group_HBS2 <- group_HBS2 %>% group_by(Sample_Name) %>% slice(which.max(time))
group_HBS2$group <- ifelse(group_HBS2$Mvalnew_bonf_mean_jmpredict <= cutoff6[1,2], "low", "high")
group_HBS2 <- group_HBS2[, c("Sample_Name", "group")]
sampledata_HBS4 <- merge(sampledata_HBS, group_HBS2, by = "Sample_Name", all.x = T)

sampledata_HBS5 <- sampledata_HBS4[,c("Sample_Name","time","status","duration","MMSE","GBA.APOE4.PHS","sex","age","educyrs","group")]
sampledata_HBS5$sex <- ifelse(sampledata_HBS5$sex == "M", "1", "2")
sampledata_HBS5$sex <- as.factor(sampledata_HBS5$sex)
sampledata_HBS5$Sample_Name <- as.factor(sampledata_HBS5$Sample_Name)
sampledata_HBS5$group <- as.factor(sampledata_HBS5$group)
sampledata_HBS5 <- sampledata_HBS5[order(sampledata_HBS5$Sample_Name, sampledata_HBS5$time),]

### build LMM for cohort
# PPMI top CpG
lme_PPMI_top1 <- lmer(moca ~ age + educyrs + duration + sex + GBA.APOE4.PHS + group1*time + (1 + time|Sample_Name), data = sampledata_PPMI2)
PPMI_top_low1 <- subset(sampledata_PPMI2, group1 == "low")
PPMI_top_high1 <- subset(sampledata_PPMI2, group1 == "high")
lme_PPMI_top_low1 <- lmer(moca ~ age + educyrs + duration + GBA.APOE4.PHS + sex + time + (1 + time|Sample_Name), data = PPMI_top_low1)
lme_PPMI_top_high1 <- lmer(moca ~ age + educyrs + duration + GBA.APOE4.PHS + sex + time + (1 + time|Sample_Name), data = PPMI_top_high1)

lme_PPMI_top2 <- lmer(moca ~ age + educyrs + duration + sex + GBA.APOE4.PHS + group2*time + (1 + time|Sample_Name), data = sampledata_PPMI2)
PPMI_top_low2 <- subset(sampledata_PPMI2, group2 == "low")
PPMI_top_high2 <- subset(sampledata_PPMI2, group2 == "high")
lme_PPMI_top_low2 <- lmer(moca ~ age + educyrs + duration + GBA.APOE4.PHS + sex + time + (1 + time|Sample_Name), data = PPMI_top_low2)
lme_PPMI_top_high2 <- lmer(moca ~ age + educyrs + duration + GBA.APOE4.PHS + sex + time + (1 + time|Sample_Name), data = PPMI_top_high2)

# PPMI integrated Mval
lme_PPMI_inter <- lmer(moca ~ age + educyrs + duration + sex + GBA.APOE4.PHS + group*time + (1 + time|Sample_Name), data = sampledata_PPMI3)
PPMI_inter_low <- subset(sampledata_PPMI3, group == "low")
PPMI_inter_high <- subset(sampledata_PPMI3, group == "high")
lme_PPMI_inter_low <- lmer(moca ~ age + educyrs + duration + GBA.APOE4.PHS + sex + time + (1 + time|Sample_Name), data = PPMI_inter_low)
lme_PPMI_inter_high <- lmer(moca ~ age + educyrs + duration + GBA.APOE4.PHS + sex + time + (1 + time|Sample_Name), data = PPMI_inter_high)


# HBS top CpG
lme_HBS_top1 <- lmer(MMSE ~ age + educyrs + duration + sex + GBA.APOE4.PHS + group1*time + (1 + time||Sample_Name), data = sampledata_HBS3)
HBS_top_low1 <- subset(sampledata_HBS3, group1 == "low")
HBS_top_high1 <- subset(sampledata_HBS3, group1 == "high")
lme_HBS_top_low1 <- lmer(MMSE ~ age + educyrs + duration + GBA.APOE4.PHS + sex + time + (1 + time||Sample_Name), data = HBS_top_low1)
lme_HBS_top_high1 <- lmer(MMSE ~ age + educyrs + duration + GBA.APOE4.PHS + sex + time + (1 + time||Sample_Name), data = HBS_top_high1)

lme_HBS_top2 <- lmer(MMSE ~ age + educyrs + duration + sex + GBA.APOE4.PHS + group2*time + (1 + time||Sample_Name), data = sampledata_HBS3)
HBS_top_low2 <- subset(sampledata_HBS3, group2 == "low")
HBS_top_high2 <- subset(sampledata_HBS3, group2 == "high")
lme_HBS_top_low2 <- lmer(MMSE ~ age + educyrs + duration + GBA.APOE4.PHS + sex + time + (1 + time||Sample_Name), data = HBS_top_low2)
lme_HBS_top_high2 <- lmer(MMSE ~ age + educyrs + duration + GBA.APOE4.PHS + sex + time + (1 + time||Sample_Name), data = HBS_top_high2)

# HBS integrated Mval
lme_HBS_inter <- lmer(MMSE ~ age + educyrs + duration + sex + GBA.APOE4.PHS + group*time + (1 + time||Sample_Name), data = sampledata_HBS5)
HBS_inter_low <- subset(sampledata_HBS5, group == "low")
HBS_inter_high <- subset(sampledata_HBS5, group == "high")
lme_HBS_inter_low <- lmer(MMSE ~ age + educyrs + duration + GBA.APOE4.PHS + sex + time + (1 + time||Sample_Name), data = HBS_inter_low)
lme_HBS_inter_high <- lmer(MMSE ~ age + educyrs + duration + GBA.APOE4.PHS + sex + time + (1 + time||Sample_Name), data = HBS_inter_high)

save(lme_PPMI_top1,lme_PPMI_top2, lme_PPMI_inter, lme_HBS_top1, lme_HBS_top2, lme_HBS_inter, lme_PPMI_top_low1,lme_PPMI_top_low2, lme_PPMI_top_high2, lme_PPMI_inter_low, lme_PPMI_inter_high, lme_HBS_top_low1, lme_HBS_top_high1, lme_HBS_top_low2, lme_HBS_top_high2, lme_HBS_inter_low, lme_HBS_inter_high, file = "PPMI_HBS_LMM.RData")

## subset P value from LMM
LMM <- c("lme_PPMI_top1", "lme_PPMI_top2", "lme_PPMI_inter", "lme_HBS_top1", "lme_HBS_top2", "lme_HBS_inter")

Pvalue <- vector("double",1)
Pvalue[1] <- summary(lme_PPMI_top1)$coefficients[9,5]
Pvalue[2] <- summary(lme_PPMI_top2)$coefficients[9,5]
Pvalue[3] <- summary(lme_PPMI_inter)$coefficients[9,5]
Pvalue[4] <- summary(lme_HBS_top1)$coefficients[9,5]
Pvalue[5] <- summary(lme_HBS_top2)$coefficients[9,5]
Pvalue[6] <- summary(lme_HBS_inter)$coefficients[9,5]

Pvalue <- as.data.frame(Pvalue)
rownames(Pvalue) <- LMM
write.csv(Pvalue, "PPMI_HBS_LMM_Pvalue.csv", row.names = T)
