setwd("D:/PPN_lab/PPN/硕士阶段/JM/JM_significant/filter/weibull_AFT_aGH_LEDD_PHS/LMM")

library(ggplot2)
library(tidyverse)

# load LMM object
Pvalue <- read.csv("PPMI_HBS_LMM_Pvalue.csv", header = T)
cutoff <- read.csv("fitJOINT_Mvalnew_cutoff6.csv", header = T)
sampledata_HBS_preLMM <- read.csv("sampledata_HBS_preLMM.csv", header = T)
sampledata_PPMI_preLMM <- read.csv("sampledata_PPMI_preLMM.csv", header = T)

# PPMI inter Mval
p3 <- 
  ggplot() + 
  geom_line(data = PPMI_pred_inter_low, aes(x = x, y = predicted), color = 'skyblue', size = 1) +          # slope
  geom_ribbon(data = PPMI_pred_inter_low, aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
              fill = "skyblue", alpha = 0.5) + # error band
  geom_line(data = PPMI_pred_inter_high, aes(x = x, y = predicted), color = 'orangered', size = 1) +          # slope
  geom_ribbon(data = PPMI_pred_inter_high, aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
              fill = "orangered", alpha = 0.5) +
  theme_classic() +
  scale_x_continuous(expand = c(0,0), limits = c(0,3)) +
  scale_y_continuous(expand = c(0,0), limits = c(13,29)) +
  xlab('Time in study (years)') +
  ylab('MoCA') +
  labs(x = "Time in study (years)",
       y = "MoCA") +
  annotate('text', x = 2.5, y = 26,
           label = Pvalue$Pvalue[3],
           size = 4, color = 'black')

p1 <- 
  ggplot(sampledata_PPMI_preLMM, aes(x = time, y = premoca, color = group ,fill = group )) +
  geom_smooth(method = "lm", se = TRUE)+ # 使用线性模型拟合数据并显示置信区间
  theme_classic()+
  scale_color_manual(values = c('orangered', 'skyblue')) +
  scale_fill_manual(values = c('orangered', 'skyblue')) +
  scale_x_continuous(expand = c(0,0), limits = c(0,3)) +
  scale_y_continuous(expand = c(0,0), limits = c(18,30)) +
  xlab('Time in study (years)') +
  ylab('Moca') +
  annotate('text', x = 2.5, y = 26,
           label = Pvalue$Pvalue[3],
           size = 4, color = 'black')

# HBS inter Mval
p4 <- 
  ggplot() + 
  geom_line(data = HBS_pred_inter_low, aes(x = x, y = predicted), color = 'skyblue', size = 1) +          # slope
  geom_ribbon(data = HBS_pred_inter_low, aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
              fill = "skyblue", alpha = 0.5) + # error band
  geom_line(data = HBS_pred_inter_high, aes(x = x, y = predicted), color = 'orangered', size = 1) +          # slope
  geom_ribbon(data = HBS_pred_inter_high, aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
              fill = "orangered", alpha = 0.5) +
  theme_classic() +
  scale_x_continuous(expand = c(0,0), limits = c(0,3)) +
  scale_y_continuous(expand = c(0,0), limits = c(24,30)) +
  xlab('Time in study (years)') +
  ylab('MoCA') +
  labs(x = "Time in study (years)",
       y = "MMSE") +
  annotate('text', x = 2.5, y = 26,
           label = Pvalue$Pvalue[6],
           size = 4, color = 'black')

p2 <- 
  ggplot(sampledata_HBS_preLMM, aes(x = time, y = preMMSE, color = group ,fill = group )) +
  geom_smooth(method = "lm", se = TRUE)+ # 使用线性模型拟合数据并显示置信区间
  theme_classic()+
  scale_color_manual(values = c('orangered', 'skyblue')) +
  scale_fill_manual(values = c('orangered', 'skyblue')) +
  scale_x_continuous(expand = c(0,0), limits = c(0,3)) +
  scale_y_continuous(expand = c(0,0), limits = c(24,30)) +
  xlab('Time in study (years)') +
  ylab('Moca') +
  annotate('text', x = 1, y = 25,
           label = Pvalue$Pvalue[6],
           size = 4, color = 'black')

#cbind the figures
library(patchwork)
library(ggpubr)
library(cowplot)
p5 <- plot_grid(p1, p2, ncol = 2, align = "vh")



