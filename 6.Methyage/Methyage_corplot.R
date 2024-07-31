setwd("D:/PPN_lab/PPN/硕士阶段/JM/JM_significant/filter/weibull_AFT_aGH_LEDD_PHS/methyage/correct")

library(ggplot2)
library(ggpubr)
library(ggpmisc)

# data input
PPMI_PCClock <- read.csv("PPMI_PCClock_DNAmAge.csv", header = T)
PPMI_PCClock2 <- PPMI_PCClock[,c(1:2,4:7,17:22)]
colnames(PPMI_PCClock2)[1] <- "ID_Position"

PPMI_sampledata <- read.csv("PPMI_sampledata_Mvalnew_methyage_correct_group.csv", header = T)
PPMI_sampledata2 <- PPMI_sampledata[,c(1:26,35)]
PPMI_sampledata3 <- merge(PPMI_sampledata2, PPMI_PCClock2, by = "ID_Position", all.x = T)

#PPMI_sampledata_BL <- subset(PPMI_sampledata3, time == 0)

HBS_PCClock <- read.csv("HBS_PCClock_DNAmAge.csv", header = T)
HBS_PCClock2 <- HBS_PCClock[,c(1:2,4:7,17:22)]
colnames(HBS_PCClock2)[1] <- "ID_Position"

HBS_sampledata <- read.csv("HBS_sampledata_Mvalnew_methyage_correct_group.csv", header = T)
HBS_sampledata2 <- HBS_sampledata[,c(1:33,39)]
HBS_sampledata3 <- merge(HBS_sampledata2, HBS_PCClock2, by = "ID_Position", all.x = T)

clockColumns = c("PCHorvath1", "PCHorvath2", "PCHannum", "PCPhenoAge", "PCGrimAge")
for (i in clockColumns){
  HBS_sampledata3[,paste0(i,"Resid2")] = HBS_sampledata3[,i] - HBS_sampledata3$Age
  PPMI_sampledata3[,paste0(i,"Resid2")] = PPMI_sampledata3[,i] - PPMI_sampledata3$Age
}


#ADNI_sampledata <- read.csv("ADNI_sampledata_Mvalnew_methyage_correct_group.csv", header = T)

### methyage ~ joint modeling score of epigenetic signature
## PPMI
# biological age
#PPMI_p1 <- ggplot(PPMI_sampledata, aes(x = Mvalnew_bonf_jmpredict, y = Age)) +
  #geom_point(size = 0.1) +
  #geom_smooth(method = "lm", color = "#8891DB", fill = "#8891DB") +
  #theme_classic() +
  #stat_cor(method = "pearson", label.x = 0, label.y = 85)

# PCHorvath1 age
#PPMI_PC1_p1 <- ggplot(PPMI_sampledata, aes(x = Mvalnew_bonf_jmpredict, y = PCHorvath1)) +
  #geom_point(size = 0.1) +
  #geom_smooth(method = "lm", color = "#E1C855", fill = "#E1C855") +
  #theme_classic() +
  #stat_cor(method = "pearson", label.x = 0, label.y = 85)

color <- c("#FCB2AF","#9BDFDF","#FAC795","#92B4C8","#BEBCDF")
#E1C855, E07B54, 51B1B7

# PCHorvath1 age acceleration
PPMI_PC1_p2 <- ggplot(PPMI_sampledata3, aes(x = Mvalnew_bonf_jmpredict, y = PCHorvath1Resid)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", color = color[1], fill = color[1]) +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0, label.y = 15) +
  xlab("Cognitive survival score") +
  ylab("Acceleration of PCHorvath1 age")

PPMI_PC1_p3 <- ggplot(PPMI_sampledata3, aes(x = Mvalnew_bonf_jmpredict, y = PCHorvath1Resid2)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", color = color[1], fill = color[1]) +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0, label.y = 15) +
  xlab("Cognitive survival score") +
  ylab("Acceleration of PCHorvath1 age")

# PCHorvath2 age
#PPMI_PC2_p1 <- ggplot(PPMI_sampledata, aes(x = Mvalnew_bonf_jmpredict, y = PCHorvath2)) +
 # geom_point(size = 0.1) +
 # geom_smooth(method = "lm", color = "#E07B54", fill = "#E07B54") +
 #  theme_classic() +
 # stat_cor(method = "pearson", label.x = 0, label.y = 85)


# PCHorvath2 age acceleration
PPMI_PC2_p2 <- ggplot(PPMI_sampledata3, aes(x = Mvalnew_bonf_jmpredict, y = PCHorvath2Resid)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", color = color[2], fill = color[2]) +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0, label.y = 16) +
  xlab("Cognitive survival score") +
  ylab("Acceleration of PCHorvath2 age")

PPMI_PC2_p3 <- ggplot(PPMI_sampledata3, aes(x = Mvalnew_bonf_jmpredict, y = PCHorvath2Resid2)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", color = color[2], fill = color[2]) +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0, label.y = 16) +
  xlab("Cognitive survival score") +
  ylab("Acceleration of PCHorvath2 age")

# PCHannum age
#PPMI_PC3_p1 <- ggplot(PPMI_sampledata, aes(x = Mvalnew_bonf_jmpredict, y = PCHannum)) +
 # geom_point(size = 0.1) +
 #geom_smooth(method = "lm", color = "#51B1B7", fill = "#51B1B7") +
  #theme_classic() +
  #stat_cor(method = "pearson", label.x = 0, label.y = 86)

# PCHannum age acceleration
PPMI_PC3_p2 <- ggplot(PPMI_sampledata3, aes(x = Mvalnew_bonf_jmpredict, y = PCHannumResid)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", color = color[3], fill = color[3]) +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0, label.y = 15) +
  xlab("Cognitive survival score") +
  ylab("Acceleration of PCHannum age")

PPMI_PC3_p3 <- ggplot(PPMI_sampledata3, aes(x = Mvalnew_bonf_jmpredict, y = PCHannumResid2)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", color = color[3], fill = color[3]) +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0, label.y = 20) +
  xlab("Cognitive survival score") +
  ylab("Acceleration of PCHannum age")

# PCPheno age acceleration
PPMI_PC4_p2 <- ggplot(PPMI_sampledata3, aes(x = Mvalnew_bonf_jmpredict, y = PCPhenoAgeResid)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", color = color[4], fill = color[4]) +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0, label.y = 15) +
  xlab("Cognitive survival score") +
  ylab("Acceleration of PCPheno age")

PPMI_PC4_p3 <- ggplot(PPMI_sampledata3, aes(x = Mvalnew_bonf_jmpredict, y = PCPhenoAgeResid2)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", color = color[4], fill = color[4]) +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0, label.y = 15) +
  xlab("Cognitive survival score") +
  ylab("Acceleration of PCPheno age")

# PCGrim age acceleration
PPMI_PC5_p2 <- ggplot(PPMI_sampledata3, aes(x = Mvalnew_bonf_jmpredict, y = PCGrimAgeResid)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", color = color[5], fill = color[5]) +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0, label.y = 10) +
  xlab("Cognitive survival score") +
  ylab("Acceleration of PCGrim age")

PPMI_PC5_p3 <- ggplot(PPMI_sampledata3, aes(x = Mvalnew_bonf_jmpredict, y = PCGrimAgeResid2)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", color = color[5], fill = color[5]) +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0, label.y = 20) +
  xlab("Cognitive survival score") +
  ylab("Acceleration of PCGrim age")

#test <- plot_grid(PPMI_PC1_p1, PPMI_PC2_p1, PPMI_PC3_p1, PPMI_p1, PPMI_PC1_p2, PPMI_PC2_p2, PPMI_PC3_p2, ncol = 4, align = "vh",
                  #rel_heights = c(1,1), rel_widths = c(1,1))
## HBS
# biological age
#HBS_p1 <- ggplot(HBS_sampledata, aes(x = Mvalnew_bonf_mean_jmpredict, y = Age)) +
 # geom_point(size = 0.1) +
  #geom_smooth(method = "lm", color = "#8891DB", fill = "#8891DB") +
  #theme_classic() +
  #stat_cor(method = "pearson", label.x = 0, label.y = 88)

# PCHorvath1 age
#HBS_PC1_p1 <- ggplot(HBS_sampledata, aes(x = Mvalnew_bonf_mean_jmpredict, y = PCHorvath1)) +
 # geom_point(size = 0.1) +
#  geom_smooth(method = "lm", color = "#E1C855", fill = "#E1C855") +
 # theme_classic() +
  #stat_cor(method = "pearson", label.x = 0, label.y = 92)

# PCHorvath1 age acceleration
HBS_PC1_p2 <- ggplot(HBS_sampledata3, aes(x = Mvalnew_bonf_mean_jmpredict, y = PCHorvath1Resid)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", color = color[1], fill = color[1]) +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0, label.y = 20) +
  xlab("Cognitive survival score") +
  ylab("Acceleration of PCHorvath1 age")

HBS_PC1_p3 <- ggplot(HBS_sampledata3, aes(x = Mvalnew_bonf_mean_jmpredict, y = PCHorvath1Resid2)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", color = color[1], fill = color[1]) +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0, label.y = 20) +
  xlab("Cognitive survival score") +
  ylab("Acceleration of PCHorvath1 age")

# PCHorvath2 age
#HBS_PC2_p1 <- ggplot(HBS_sampledata, aes(x = Mvalnew_bonf_mean_jmpredict, y = PCHorvath2)) +
 # geom_point(size = 0.1) +
  #geom_smooth(method = "lm", color = "#E07B54", fill = "#E07B54") +
  #theme_classic() +
  #stat_cor(method = "pearson", label.x = 0, label.y = 85)

# PCHorvath2 age acceleration
HBS_PC2_p2 <- ggplot(HBS_sampledata3, aes(x = Mvalnew_bonf_mean_jmpredict, y = PCHorvath2Resid)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", color = color[2], fill = color[2]) +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0, label.y = 20) +
  xlab("Cognitive survival score") +
  ylab("Acceleration of PCHorvath2 age")

HBS_PC2_p3 <- ggplot(HBS_sampledata3, aes(x = Mvalnew_bonf_mean_jmpredict, y = PCHorvath2Resid2)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", color = color[2], fill = color[2]) +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0, label.y = 20) +
  xlab("Cognitive survival score") +
  ylab("Acceleration of PCHorvath2 age")

# PCHannum age
#HBS_PC3_p1 <- ggplot(HBS_sampledata, aes(x = Mvalnew_bonf_mean_jmpredict, y = PCHannum)) +
 # geom_point(size = 0.1) +
  #geom_smooth(method = "lm", color = "#51B1B7", fill = "#51B1B7") +
  #theme_classic() +
  #stat_cor(method = "pearson", label.x = 0, label.y = 92)

# PCHannum age acceleration
HBS_PC3_p2 <- ggplot(HBS_sampledata3, aes(x = Mvalnew_bonf_mean_jmpredict, y = PCHannumResid)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", color = color[3], fill = color[3]) +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0, label.y = 20) +
  xlab("Cognitive survival score") +
  ylab("Acceleration of PCHannum age")

HBS_PC3_p3 <- ggplot(HBS_sampledata3, aes(x = Mvalnew_bonf_mean_jmpredict, y = PCHannumResid2)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", color = color[3], fill = color[3]) +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0, label.y = 20) +
  xlab("Cognitive survival score") +
  ylab("Acceleration of PCHannum age")

# PCPheno age acceleration
HBS_PC4_p2 <- ggplot(HBS_sampledata3, aes(x = Mvalnew_bonf_mean_jmpredict, y = PCPhenoAgeResid)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", color = color[4], fill = color[4]) +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0, label.y = 23) +
  xlab("Cognitive survival score") +
  ylab("Acceleration of PCPheno age")

HBS_PC4_p3 <- ggplot(HBS_sampledata3, aes(x = Mvalnew_bonf_mean_jmpredict, y = PCPhenoAgeResid2)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", color = color[4], fill = color[4]) +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0, label.y = 23) +
  xlab("Cognitive survival score") +
  ylab("Acceleration of PCPheno age")

# PCGrim age acceleration
HBS_PC5_p2 <- ggplot(HBS_sampledata3, aes(x = Mvalnew_bonf_mean_jmpredict, y = PCGrimAgeResid)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", color = color[5], fill = color[5]) +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0, label.y = 10) +
  xlab("Cognitive survival score") +
  ylab("Acceleration of PCGrim age")

HBS_PC5_p3 <- ggplot(HBS_sampledata3, aes(x = Mvalnew_bonf_mean_jmpredict, y = PCGrimAgeResid2)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", color = color[5], fill = color[5]) +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0, label.y = 20) +
  xlab("Cognitive survival score") +
  ylab("Acceleration of PCGrim age")


## ADNI
# biological age
#ADNI_p1 <- ggplot(ADNI_sampledata, aes(x = Mvalnew_bonf_jmpredict, y = Age)) +
 # geom_point(size = 0.1) +
  #geom_smooth(method = "lm", color = "#8891DB", fill = "#8891DB") +
  #theme_classic() +
  #stat_cor(method = "pearson", label.x = 0, label.y = 92)

# PCHorvath1 age
#ADNI_PC1_p1 <- ggplot(ADNI_sampledata, aes(x = Mvalnew_bonf_jmpredict, y = PCHorvath1)) +
 # geom_point(size = 0.1) +
  #geom_smooth(method = "lm", color = "#E1C855", fill = "#E1C855") +
  #theme_classic() +
  #stat_cor(method = "pearson", label.x = 0, label.y = 92)

# PCHorvath1 age acceleration
#ADNI_PC1_p2 <- ggplot(ADNI_sampledata, aes(x = Mvalnew_bonf_jmpredict, y = PCHorvath1Resid)) +
 # geom_point(size = 0.1) +
  #geom_smooth(method = "lm", color = "#E1C855", fill = "#E1C855") +
  #theme_classic() +
  #stat_cor(method = "pearson", label.x = 0, label.y = 20)

# PCHorvath2 age
#ADNI_PC2_p1 <- ggplot(ADNI_sampledata, aes(x = Mvalnew_bonf_jmpredict, y = PCHorvath2)) +
#geom_point(size = 0.1) +
#geom_smooth(method = "lm", color = "#E07B54", fill = "#E07B54") +
#theme_classic() +
#stat_cor(method = "pearson", label.x = 0, label.y = 92)

# PCHorvath2 age acceleration
#ADNI_PC2_p2 <- ggplot(ADNI_sampledata, aes(x = Mvalnew_bonf_jmpredict, y = PCHorvath2Resid)) +
#geom_point(size = 0.1) +
#geom_smooth(method = "lm", color = "#E07B54", fill = "#E07B54") +
#theme_classic() +
#stat_cor(method = "pearson", label.x = 0, label.y = 30)

# PCHannum age
#ADNI_PC3_p1 <- ggplot(ADNI_sampledata, aes(x = Mvalnew_bonf_jmpredict, y = PCHannum)) +
# geom_point(size = 0.1) +
#geom_smooth(method = "lm", color = "#51B1B7", fill = "#51B1B7") +
#theme_classic() +
#stat_cor(method = "pearson", label.x = 0, label.y = 100)

# PCHannum age acceleration
#ADNI_PC3_p2 <- ggplot(ADNI_sampledata, aes(x = Mvalnew_bonf_jmpredict, y = PCHannumResid)) +
#geom_point(size = 0.1) +
#geom_smooth(method = "lm", color = "#51B1B7", fill = "#51B1B7") +
# theme_classic() +
#stat_cor(method = "pearson", label.x = 0, label.y = 28)

# cbind the figures
library(patchwork)
library(cowplot)

# include biological age and acceleration of different PC clock ages from PPMI and HBS cohort
corplot1 <- plot_grid(PPMI_PC1_p2, PPMI_PC2_p2, PPMI_PC3_p2, PPMI_PC4_p2,PPMI_PC5_p2,HBS_PC1_p2, HBS_PC2_p2, HBS_PC3_p2,HBS_PC4_p2,HBS_PC5_p2,
                     ncol = 5, align = "vh", rel_heights = c(1,1), rel_widths = c(1,1))

corplot2 <- plot_grid(PPMI_PC1_p3, PPMI_PC2_p3, PPMI_PC3_p3, PPMI_PC4_p3,PPMI_PC5_p3,HBS_PC1_p3, HBS_PC2_p3, HBS_PC3_p3,HBS_PC4_p3,HBS_PC5_p3,
                      ncol = 5, align = "vh", rel_heights = c(1,1), rel_widths = c(1,1))
# include different PC clock ages from PPMI and HBS
#corplot2 <- plot_grid(PPMI_PC1_p1, PPMI_PC2_p1, PPMI_PC3_p1, HBS_PC1_p1, HBS_PC2_p1, HBS_PC3_p1, ncol = 3, align = "vh",
# rel_heights = c(1,1), rel_widths = c(1,1))
# include different PC clock ages, acceleration and biological age from ADNI
#corplot3 <- plot_grid(ADNI_PC1_p1, ADNI_PC2_p1, ADNI_PC3_p1, ADNI_PC1_p2, ADNI_PC2_p2, ADNI_PC3_p2, ncol = 3,
#align = "vh", rel_heights = c(1,1), rel_widths = c(1,1))
# include all ages from all cohorts
#corplot4 <- plot_grid(PPMI_p1, PPMI_PC1_p1, PPMI_PC2_p1, PPMI_PC3_p1,
# HBS_p1, HBS_PC1_p1, HBS_PC2_p1, HBS_PC3_p1,
#ADNI_p1, ADNI_PC1_p1, ADNI_PC2_p1, ADNI_PC3_p1,
# ncol = 4, align = "vh", rel_heights = c(1,1), rel_widths = c(1,1))
#colplot5 <- plot_grid(PPMI_p1, PPMI_PC1_p2, PPMI_PC2_p2, PPMI_PC3_p2,
# HBS_p1, HBS_PC1_p2, HBS_PC2_p2, HBS_PC3_p2,
#ADNI_p1, ADNI_PC1_p2, ADNI_PC2_p2, ADNI_PC3_p2,
#ncol = 4, align = "vh", rel_heights = c(1,1), rel_widths = c(1,1))

## methyage ~ biological age
formula <- y ~ x
# Horvath1
PPMI_p1 <- ggplot(PPMI_sampledata, aes(x = Age, y = PCHorvath1)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", color = "#E1C855", fill = "#E1C855") +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 40, label.y = 78)
  #stat_poly_eq(aes(label = ..eq.label..), formula = formula,parse = TRUE, geom = "text", label.x = 50, label.y = 76, hjust = 0)

HBS_p1 <- ggplot(HBS_sampledata, aes(x = Age, y = PCHorvath1)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", color = "#E1C855", fill = "#E1C855") +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 60, label.y = 85)

#Horvath2
PPMI_p2 <- ggplot(PPMI_sampledata, aes(x = Age, y = PCHorvath2)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", color = "#E07B54", fill = "#E07B54") +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 40, label.y = 78) 

HBS_p2 <- ggplot(HBS_sampledata, aes(x = Age, y = PCHorvath2)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", color = "#E07B54", fill = "#E07B54") +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 60, label.y = 80)

#Hannum
PPMI_p3 <- ggplot(PPMI_sampledata, aes(x = Age, y = PCHannum)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", color = "#51B1B7", fill = "#51B1B7") +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 40, label.y = 78)

HBS_p3 <- ggplot(HBS_sampledata, aes(x = Age, y = PCHannum)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", color = "#51B1B7", fill = "#51B1B7") +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 60, label.y = 85)



p4 <- plot_grid(PPMI_p1, PPMI_p2, PPMI_p3, HBS_p1, HBS_p2, HBS_p3, ncol = 3, align = "vh",
                rel_heights = c(1,1), rel_widths = c(1,1))
