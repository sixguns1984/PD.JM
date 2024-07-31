
setwd("D:/PPN_lab/PPN/JM/JM_significant/filter/Volcano")
piecewise_PH_aGH <- read.csv("fitJOINT_piecewise_PH_aGH_scale_merge.csv",header = T)
spline_PH_aGH <- read.csv("fitJOINT_spline_PH_aGH_scale_FDR.csv",header = T)
weibull_AFT_aGH <- read.csv("fitJOINT_weibull_AFT_aGH_scale_FDR.csv",header = T)
weibull_PH_aGH <- read.csv("fitJOINT_weibull_PH_aGH_scale_merge.csv",header = T)

#org.Hs.eg.db annotation
CpG_annotation_peak_Anno <- data.table::fread("D:/PPN_lab/PPN/JM/methylation_annotation/CpG_annotation_peak_Anno.txt",sep = '\t')
CpG_annotation <- CpG_annotation_peak_Anno[,c(1,6,7,13,14)]
CpG_annotation$geneId <- substring(CpG_annotation$geneId,1,15)
names(CpG_annotation)[4] <- "ENSEMBL"

#genecode annotation, more intact
library(rtracklayer)
genecode <- rtracklayer::import('D:/PPN_lab/PPN/JM/methylation_annotation/gencode.v42.basic.annotation.gtf')
genecode <- as.data.frame(genecode)
genecode_annotaiton <- genecode[,c(10:12)]
genecode_annotaiton$gene_id <- substring(genecode_annotaiton$gene_id,1,15)
names(genecode_annotaiton)[1] <- "ENSEMBL"
genecode_annotaiton <- genecode_annotaiton[!duplicated(genecode_annotaiton$ENSEMBL),]
CpG_genecode_annotation <- merge(CpG_annotation, genecode_annotaiton, by = 'ENSEMBL', all.x = T)
names(CpG_genecode_annotation)[7] <- "SYMBOL"

piecewise_PH_aGH_annotation_genecode <- merge(piecewise_PH_aGH, CpG_genecode_annotation, by = "Name", all.x = T)
piecewise_PH_aGH_annotation <- piecewise_PH_aGH_annotation_genecode[order(piecewise_PH_aGH_annotation_genecode$pvalue_FDR.y, decreasing = F),]

weibull_PH_aGH_annotation_genecode <- merge(weibull_PH_aGH, CpG_genecode_annotation, by = 'Name', all.x = T)
weibull_PH_aGH_annotation <- weibull_PH_aGH_annotation_genecode[order(weibull_PH_aGH_annotation_genecode$pvalue_FDR, decreasing = F),]

spline_PH_aGH_annotation_genecode <- merge(spline_PH_aGH, CpG_genecode_annotation, by = 'Name', all.x = T)
spline_PH_aGH_annotation <- spline_PH_aGH_annotation_genecode[order(spline_PH_aGH_annotation_genecode$pvalue_FDR, decreasing = F),]

weibull_AFT_aGH_annotation_genecode <- merge(weibull_AFT_aGH, CpG_genecode_annotation, by = 'Name', all.x = T)
weibull_AFT_aGH_annotation <- weibull_AFT_aGH_annotation_genecode[order(weibull_AFT_aGH_annotation_genecode$pvalue_FDR, decreasing = F),]

#volcano plot (ggplot2)
library(ggplot2)
library(ggrepel)
data1 <- piecewise_PH_aGH_annotation
data2 <- weibull_PH_aGH_annotation
data3 <- spline_PH_aGH_annotation
data4 <- weibull_AFT_aGH_annotation
data1$changed <- factor(ifelse(data1$pvalue_bonferroni.y < 0.05,'yes','no')) #Screen for CpGs with significant changes
data2$changed <- factor(ifelse(data2$pvalue_bonferroni < 0.05,'yes','no'))
data3$changed <- factor(ifelse(data3$pvalue_bonferroni < 0.05,'Bonferroni < 0.05','Bonferroni ≥ 0.05'))
data3$changed <- factor(data3$changed, levels=c('Bonferroni < 0.05','Bonferroni ≥ 0.05'))
data4$changed <- factor(ifelse(data4$pvalue_bonferroni < 0.05,'yes','no'))

rm(CpG_annotation_peak_Anno, genecode, test)

p1 <- ggplot(data1, aes(abeta.y, -log10(pvalue.y), color = factor(changed))) +  
        geom_point(size = 2) +
        labs(x = NULL, y = NULL, title = "Piecewise-PH", hjust=0.5) +
        theme_classic(base_size = 12) +
        scale_color_manual(values = c('#9ca8b8','#c4403b')) +
        theme(legend.position = "none", 
          plot.title = element_text(hjust=0.5, face="bold")) +
        geom_text_repel(data = data1[data1$pvalue_bonferroni.y < 0.05,], aes(label = SYMBOL), color = "#c4403b", size = 3, max.overlaps = 50, #add text tags to significant CpGs
                  box.padding = unit(0.5, "lines"), 
                  point.padding = NA, 
                  segment.colour = "black") +
        theme(axis.text= element_text(colour = "black"),
            panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
        expand_limits(x=c(-20,20), y=c(0, 10))

p2 <- ggplot(data2, aes(abeta.x, -log10(pvalue.x), color = factor(changed))) +  
        geom_point(size = 2) +
        labs(x = NULL, y = NULL, title = "Weibull-PH", hjust=0.5) +
        theme_classic(base_size = 12) +
        scale_color_manual(values = c('#9ca8b8','#c4403b')) +
        theme(legend.position = "none", 
          plot.title = element_text(hjust=0.5, face="bold")) +
        geom_text_repel(data = data2[data2$pvalue_bonferroni < 0.05,], aes(label = SYMBOL), color = "#c4403b", size = 3, max.overlaps = 50, #add text tags to significant CpGs
                  box.padding = unit(0.5, "lines"), 
                  point.padding = NA, 
                  segment.colour = "black") +
        theme(axis.text= element_text(colour = "black"),
            panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
        coord_cartesian(xlim = c(-12, 12), ylim = c(0, 8))

p3 <- ggplot(data3, aes(abeta, -log10(pvalue), color = factor(changed))) +  
        geom_point(size = 2) +
        labs(x = NULL, y = NULL, title = "Spline-PH", hjust=0.5) +
        theme_classic(base_size = 12) +
        scale_color_manual(values = c('#c4403b', '#9ca8b8')) +
        theme(legend.title = element_blank(),
              legend.position = c(0.85,0.95), 
          plot.title = element_text(hjust=0.5, face="bold")) +
        geom_text_repel(data = data3[data3$pvalue_bonferroni < 0.05,], aes(label = SYMBOL), color = "#c4403b", size = 3, max.overlaps = 50, #add text tags to significant CpGs
                  box.padding = unit(0.5, "lines"), 
                  point.padding = NA, 
                  segment.colour = "black") +
        theme(axis.text= element_text(colour = "black"),
            panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
        coord_cartesian(xlim = c(-20, 20), ylim = c(0, 12))
  
p4 <- ggplot(data4, aes(abeta, -log10(pvalue), color = factor(changed))) +  
        geom_point(size = 2) +
        labs(x = NULL, y = NULL, title = "Weibull-AFT", hjust=0.5) +
        theme_classic(base_size = 12) +
        scale_color_manual(values = c('#9ca8b8','#c4403b')) +
        theme(legend.position = "none", 
          plot.title = element_text(hjust=0.5, face="bold")) +
        geom_text_repel(data = data4[data4$pvalue_bonferroni < 0.05,], aes(label = SYMBOL), color = "#c4403b", size = 3, max.overlaps = 50, #add text tags to significant CpGs
                  box.padding = unit(0.5, "lines"), 
                  point.padding = NA, 
                  segment.colour = "black") +
        theme(axis.text= element_text(colour = "black"),
            panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
        coord_cartesian(xlim = c(-2.5, 2.5), ylim = c(0, 8))



#cbind the figures
library(patchwork)
library(ggpubr)
library(cowplot)

p5 <- plot_grid(p1, p3, p2, p4, ncol = 2, align = "vh")
annotate_figure(p5, left = text_grob(expression(-Log[10]*" (P-value)"), rot = 90, vjust = 0.5, size = 14),
                bottom = text_grob("Effect size", hjust = 0.5, size = 14))
