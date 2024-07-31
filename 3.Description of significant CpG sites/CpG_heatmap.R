setwd("D:/PPN_lab/PPN/硕士阶段/JM/JM_significant/filter/weibull_AFT_aGH_LEDD_PHS/methylation")

## data prepared
library(dplyr)
fitJOINT_annotation <- read.csv("fitJOINT_weibull_AFT_annotation.csv", header = T)
fitJOINT005_annotation <- subset(fitJOINT_annotation, pvalue_bonferroni < 0.05)

sampledata <- read.csv("sampledata_bonf.csv", header = T)
sampledata$Sample_Group <- as.factor(sampledata$Sample_Group)
sampledata$sex <- as.factor(sampledata$sex)
sampledata$status <- as.factor(sampledata$status)


## multiple time
sampledata <- sampledata[order(sampledata$time,sampledata$Sample_Name),]
rownames(sampledata) <- paste(sampledata$Sample_Name, sampledata$time, sep = "_")
test4 <- sampledata[,-c(1:6)]
test5 <- t(test4)
test5 <- as.data.frame(test5)

annotation_col2 <- as.data.frame(sampledata[,c("time","status")])
rownames(annotation_col2) <- rownames(sampledata)
names(annotation_col2) <- c("time","status")
annotation_col2$status <- as.factor(annotation_col2$status)
annotation_col2$time <- as.factor(annotation_col2$time)

annotation_row <- fitJOINT005_annotation[,c(1,10)]
annotation_row2 <- annotation_row %>% mutate(Annotation = ifelse(grepl("Promoter", annotation_row$annotation), "Promoter", 
                                                                 ifelse(grepl("Intron", annotation_row$annotation), "Intron", 
                                                                        ifelse(grepl("Exon", annotation_row$annotation), "Exon", 
                                                                               ifelse(grepl("3' UTR", annotation_row$annotation), "3' UTR",
                                                                                      ifelse(grepl("5' UTR", annotation_row$annotation), "5' UTR", "Distal Intergenic"))))))
annotation_row3 <- as.data.frame(annotation_row2[, -c(1:2)])
rownames(annotation_row3) <- annotation_row2$Name
names(annotation_row3) <- "annotation"
annotation_row3$annotation <- as.factor(annotation_row3$annotation)

## Time point clustering is performed on columns
time_point_order <- c("0","1","2","3")
time_points <- factor(annotation_col2$time, levels = time_point_order)

## Computes the row and column distance matrix
row_dist <- as.dist(1 - cor(test5))

## Computes row and column hierarchical clustering
row_clusters <- hclust(row_dist)


pheatmap(test5, annotation_row = annotation_row3,
         cluster_rows = TRUE, cluster_cols = FALSE,
         gaps_col = c(171, 324, 445), show_colnames = FALSE, 
         angle_row = 0, annotation_names_row = FALSE,
         annotation_col = ,
         clustering_distance_cols = list(
           BL = "euclidean",
           V04 = "euclidean",
           V06 = "euclidean",
           V08 = "euclidean")) 
         
