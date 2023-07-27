library('rtracklayer')
library('dplyr')
library("reshape2")
library("ggplot2")

setwd("X:/Shamitha/tinyrna/START_HERE/results/tinyrna_2023-07-17_20-17-19_Cel_2021_08_GA_26G/tiny-deseq")
deseq_ga<- read.csv("tinyrna_2023-07-17_20-17-19_cond1_Cel_Del_GA_cond2_Cel_WT_GA_deseq_table.csv")
ggplot(deseq_ga, aes(y = log2FoldChange, x = log2(baseMean))) +
  geom_point(color = "grey", size = 1) +
  theme_minimal() +
  labs(y = "Log2 Fold Change", x = "Log2 Mean of Counts", title = "Differentially targeted genes in C. elegans Gravid Adults (p<0.01)") +
  geom_point(data = subset(deseq_ga, deseq_ga$padj<0.01), color = "#F07167", size = 1)+
  geom_point(data = subset(deseq_ga, deseq_ga$Classifier %in% c("Genes 26G ", "Pseudogene 26G", "Pseudogene 26G") & deseq_ga$log2FoldChange< -5), color = "#009975", size = 1)


setwd("X:/Shamitha/tinyrna/START_HERE/results/tinyrna_2023-07-17_21-29-34_Cel_2021_08_L4_26G/tiny-deseq")
deseq_l4<- read.csv("tinyrna_2023-07-17_21-29-34_cond1_Cel_Del_L4_cond2_Cel_WT_L4_deseq_table.csv")
ggplot(deseq_l4, aes(y = log2FoldChange, x = log2(baseMean))) +
  geom_point(color = "grey", size = 1) +
  theme_minimal() +
  labs(y = "Log2 Fold Change", x = "Log2 Mean of Counts", title = "Differentially targeted genes in C. elegans L4 (p<0.01)") +
  geom_point(data = subset(deseq_l4, deseq_l4$padj<0.01), color = "#F07167", size = 1)+
  geom_point(data = subset(deseq_l4, deseq_l4$Classifier %in% c("Genes 26G ", "Pseudogene 26G", "Pseudogene 26G") & deseq_l4$log2FoldChange< -5), color = "#009975", size = 1)


setwd("X:/Shamitha/tinyrna/START_HERE/tinyrna_2023-07-18_15-19-45_cbr_26G_GA/tiny-deseq")
deseq_ga<- read.csv("tinyrna_2023-07-18_15-19-45_cond1_gtsf-1(xf345)_GA_cond2_WT_GA_deseq_table.csv")
ggplot(deseq_ga, aes(y = log2FoldChange, x = log2(baseMean))) +
  geom_point(color = "grey", size = 1) +
  theme_minimal() +
  labs(y = "Log2 Fold Change", x = "Log2 Mean of Counts", title = "Differentially targeted genes in C. briggsae Gravid Adults (p<0.01)") +
  geom_point(data = subset(deseq_ga, deseq_ga$padj<0.01), color = "#F07167", size = 1)+
  geom_point(data = subset(deseq_ga, deseq_ga$Classifier %in% c("Genes 26G ", "Pseudogene 26G", "Pseudogene 26G") & deseq_ga$log2FoldChange< -5), color = "#009975", size = 1)

setwd("X:/Shamitha/tinyrna/START_HERE/tinyrna_2023-07-18_14-02-56_cbr_26G_L4/tiny-deseq")
deseq_l4<- read.csv("tinyrna_2023-07-18_14-02-56_cond1_gtsf-1(xf345)_L4_cond2_WT_L4_deseq_table.csv")
ggplot(deseq_l4, aes(y = log2FoldChange, x = log2(baseMean))) +
  geom_point(color = "grey", size = 1) +
  theme_minimal() +
  labs(y = "Log2 Fold Change", x = "Log2 Mean of Counts", title = "Differentially targeted genes in C. briggsae L4 (p<0.01)") +
  geom_point(data = subset(deseq_l4, deseq_l4$padj<0.01), color = "#F07167", size = 1)+
  geom_point(data = subset(deseq_l4, deseq_l4$Classifier %in% c("Genes 26G ", "Pseudogene 26G", "Pseudogene 26G") & deseq_l4$log2FoldChange< -5), color = "#009975", size = 1)




