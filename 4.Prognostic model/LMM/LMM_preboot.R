setwd("D:/PPN_lab/PPN/硕士阶段/JM/JM_significant/filter/weibull_AFT_aGH_LEDD_PHS/LMM")

library(ggplot2)
library(tidyverse)
library(lme4)
library(lmerTest)

Pvalue <- read.csv("PPMI_HBS_LMM_Pvalue.csv", header = T)
cutoff <- read.csv("fitJOINT_Mvalnew_cutoff6.csv", header = T)
sampledata_PPMI <- read.csv("sampledata_PPMI_Mvalnew.csv", header = T)
sampledata_HBS <- read.csv("sampledata_HBS_preLMM.csv", header = T)

## grouping
group_PPMI = sampledata_PPMI[order(sampledata_PPMI$Sample_Name, sampledata_PPMI$time),]
group_PPMI <- group_PPMI %>% group_by(Sample_Name) %>% slice(which.max(time))
group_PPMI$group <- ifelse(group_PPMI$Mvalnew_bonf_jmpredict <= cutoff[1,2], "low", "high")
group_PPMI <- group_PPMI[, c("Sample_Name", "group")]
sampledata_PPMI2 <- merge(sampledata_PPMI, group_PPMI, by = "Sample_Name", all.x = T)

## data process
sampledata_PPMI2$age <- 10^(sampledata_PPMI2$age)
sampledata_PPMI2$sex <- as.factor(sampledata_PPMI2$sex)
sampledata_PPMI2$Sample_Name <- as.character(sampledata_PPMI2$Sample_Name)
sampledata_PPMI2$group <- as.factor(sampledata_PPMI2$group)
sampledata_PPMI2 <- sampledata_PPMI2[order(sampledata_PPMI2$Sample.Event, decreasing = F),]


sampledata_HBS$group <- as.factor(sampledata_HBS$group)
sampledata_HBS$sex <- as.factor(sampledata_HBS$sex)


## LMM
# PPMI
lme_PPMI <- lmer(moca ~ age + educyrs + duration + sex + GBA.APOE4.PHS + group*time + (1 + time||Sample_Name), data = sampledata_PPMI2)
summary(lme_PPMI)

library(bootpredictlme4)
options(bootnsim = 1000)
newdata_PPMI <- predict(lme_PPMI, newdata = sampledata_PPMI2, re.form = NA, se.fit = TRUE)
sampledata_PPMI2$predict = newdata_PPMI$fit
#sampledata_PPMI2$selo <- sampledata_PPMI2$predict - newdata_PPMI$se.fit
#sampledata_PPMI2$sehi <- sampledata_PPMI2$predict + newdata_PPMI$se.fit


# HBS
lme_HBS <- lmer(MMSE ~ age + educyrs + duration + sex + GBA.APOE4.PHS + group*time + (1 + time||Sample_Name), data = sampledata_HBS)
summary(lme_HBS)

newdata_HBS <- predict(lme_HBS, newdata=sampledata_HBS, re.form=NA, se.fit=TRUE, bootnsim = 1000)
sampledata_HBS$predict <- newdata_HBS$fit
#sampledata_HBS$selo <- sampledata_HBS$predict - newdata_HBS$se.fit
#sampledata_HBS$sehi <- sampledata_HBS$predict + newdata_HBS$se.fit

## plot
# PPMI
P1 <- ggplot(data=sampledata_PPMI2, aes(x=time, y=predict, color=group, fill = group)) +
  geom_smooth(method = "lm", se = TRUE) +
  #geom_ribbon(data=sampledata_PPMI2, 
              #aes(x=time, y=predict, ymin=selo, ymax=sehi, fill=group), alpha = 0.5)+
  theme_classic()+
  scale_color_manual(values = c('orangered', 'skyblue')) +
  scale_fill_manual(values = c('orangered', 'skyblue')) +
  scale_x_continuous(expand = c(0,0), limits = c(0,3)) +
  scale_y_continuous(expand = c(0,0), limits = c(16,30), n.breaks = 6) +
  xlab('Time in study (years)') +
  ylab('MoCA') +
  annotate('text', x = 2.5, y = 26,
           label = Pvalue$Pvalue[3],
           size = 4, color = 'black')


# HBS
P2 <- ggplot(data=sampledata_HBS, aes(x=time, y=predict, color=group, fill = group)) +
  geom_smooth(method = "lm", se = TRUE) +
  #geom_ribbon(data=sampledata_HBS, 
  #aes(x=time, y=predict, ymin=selo, ymax=sehi, fill=group), alpha = 0.5)+
  theme_classic()+
  scale_color_manual(values = c('orangered', 'skyblue')) +
  scale_fill_manual(values = c('orangered', 'skyblue')) +
  scale_x_continuous(expand = c(0,0), limits = c(0,3)) +
  scale_y_continuous(expand = c(0,0), limits = c(24,30)) +
  xlab('Time in study (years)') +
  ylab('MMSE') +
  annotate('text', x = 2, y = 26,
           label = Pvalue$Pvalue[6],
           size = 4, color = 'black')

#cbind the figures
library(patchwork)
library(ggpubr)
library(cowplot)
P3 <- plot_grid(P1, P2, ncol = 2, align = "vh")
