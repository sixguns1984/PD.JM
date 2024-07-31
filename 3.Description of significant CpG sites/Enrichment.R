setwd("D:/PPN_lab/PPN/硕士阶段/JM/JM_significant/filter/weibull_AFT_aGH_LEDD_PHS/enrichment")

fitJOINT_ann <- read.csv("fitJOINT_weibull_AFT_annotation.csv", header = T)
CpG_005 <- subset(fitJOINT_ann, pvalue_bonferroni < 0.05)

## GO enrichment
library(methylGSA)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
pvalue <- CpG_005$pvalue
names(pvalue) <- CpG_005$Name
enrichGO <- methylglm(cpg.pval = pvalue, array.type = "EPIC", minsize = 100, 
                 maxsize = 1000, GS.type = "GO", group = "all", GS.idtype = "SYMBOL")

barplot(enrichGO, num = 15, colorby = "pvalue")
write.csv(enrichGO, "enrichGO_GSA_glm.csv", row.names = F)

# KEGG enrichment
enrichKEGG <- methylglm(cpg.pval = pvalue, array.type = "EPIC", minsize = 100,
                        maxsize = 1000, GS.type = "KEGG", group = "all", GS.idtype = "ENSEMBL")

barplot(enrichKEGG, colorby = "pvalue")
write.csv(enrichKEGG, "enrichKEGG_GSA_glm.csv", row.names = F)

## GSEA
library(ChAMP)
sig_Mval <- read.csv("FDR_Mval_unscale.csv", header = T)
enrichGSEA <- champ.GSEA(beta = Mval2, DMP = NULL, DMR = NULL, 
                         CpGlist = sig$Name, method = "fisher",
                         arraytype = "EPIC", Rplot = FALSE, adjPval = 0.05)
## GO enrichment
library(missMethyl)
enrichGO2 <- gometh(sig.cpg = CpG_005$Name, collection = "GO",
                    array.type = "EPIC", sig.genes = TRUE)

write.csv(enrichGO2, "enrichGO_miss_gometh.csv", row.names = F)

## KEGG enrichment
enrichKEGG2 <- gometh(sig.cpg = CpG_005$Name, collection = "KEGG",
                      array.type = "EPIC", sig.genes = TRUE)
write.csv(enrichKEGG2, "enrichKEGG_miss_gometh.csv", row.names = F)

### visualization
# GO dotplot
library(ggplot2)
enrichGO2 <- read.csv("enrichGO_miss_gometh.csv", header = T)
enrichGO2 <- enrichGO2[order(enrichGO2$P.DE, decreasing = F), ]
GO_BP <- subset(enrichGO2, ONTOLOGY == "BP" & P.DE < 0.005)
GO_MF <- subset(enrichGO2, ONTOLOGY == "MF" & P.DE < 0.01)
GO_CC <- subset(enrichGO2, ONTOLOGY == "CC" & P.DE < 0.05)
enrichGO3 <- rbind(GO_BP, GO_MF, GO_CC)
enrichGO3$DE <- as.factor(enrichGO3$DE)
names(enrichGO3)[4] <- "Count"
names(enrichGO3)[5] <- "Pvalue"
enrichGO3$Count <- as.numeric(enrichGO3$Count)
enrichGO3$CpGRatio <- enrichGO3$Count / enrichGO3$N

ggplot(enrichGO3, aes(x = CpGRatio, y = TERM)) +
  geom_point(aes(size = Count, color = Pvalue)) +
  labs(y = "GO Term", title = "GO Enrichment", 
       color = expression(Pvalue, size = "Count") ) +
  scale_color_gradient(low = "red", high = "blue") +
  facet_grid(ONTOLOGY~., scale = 'free_y', space = "free_y") +
  theme_bw()

# KEGG
library(enrichplot)
enrichKEGG2 <- read.csv("enrichKEGG_miss_gometh.csv", header = T)
enrichKEGG3 <- subset(enrichKEGG2, P.DE < 0.05)
cnetplot(enrichKEGG3)
