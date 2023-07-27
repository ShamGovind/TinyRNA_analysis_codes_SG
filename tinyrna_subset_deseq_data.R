library('rtracklayer')
library('dplyr')
library("reshape2")
library("ggplot2")

setwd("X:/Shamitha/tinyrna/START_HERE/results/tinyrna_2023-07-17_20-17-19_Cel_2021_08_GA_26G/tiny-deseq")
deseq_ga<- read.csv("tinyrna_2023-07-17_20-17-19_cond1_Cel_Del_GA_cond2_Cel_WT_GA_deseq_table.csv")
df<-subset(deseq_ga,deseq_ga$log2FoldChange< (-5) & deseq_ga$padj<0.01)
unique(df$Classifier)
write.csv(df, "X:/Shamitha/tinyrna/START_HERE/analysis/targets_26G_cel_GA_hc.csv")


setwd("X:/Shamitha/tinyrna/START_HERE/results/tinyrna_2023-07-17_21-29-34_Cel_2021_08_L4_26G/tiny-deseq")
deseq_l4<- read.csv("tinyrna_2023-07-17_21-29-34_cond1_Cel_Del_L4_cond2_Cel_WT_L4_deseq_table.csv")
df<-subset(deseq_l4,deseq_l4$log2FoldChange< (-5) & deseq_l4$padj<0.01)
unique(df$Classifier)
write.csv(df, "X:/Shamitha/tinyrna/START_HERE/analysis/gtsf1_targets_cel_L4.csv")




setwd("X:/Shamitha/tinyrna/START_HERE/tinyrna_2023-07-18_15-19-45_cbr_26G_GA/tiny-deseq")
deseq_GA<- read.csv("tinyrna_2023-07-18_15-19-45_cond1_gtsf-1(xf345)_GA_cond2_WT_GA_deseq_table.csv")
df<-subset(deseq_GA,deseq_GA$log2FoldChange< (-5) & deseq_GA$padj<0.01)
unique(df$Classifier)
write.csv(df, "X:/Shamitha/tinyrna/START_HERE/analysis/gtsf1_targets_cbg_GA.csv")

setwd("X:/Shamitha/tinyrna/START_HERE/tinyrna_2023-07-18_14-02-56_cbr_26G_L4/tiny-deseq")
deseq_l4<- read.csv("tinyrna_2023-07-18_14-02-56_cond1_gtsf-1(xf345)_L4_cond2_WT_L4_deseq_table.csv")
df<-subset(deseq_l4,deseq_l4$log2FoldChange< (-5) & deseq_l4$padj<0.01)
unique(df$Classifier)
write.csv(df, "X:/Shamitha/tinyrna/START_HERE/analysis/gtsf1_targets_cbg_L4.csv")




setwd( "X:/Shamitha/tinyrna/START_HERE/tinyrna_2023-07-21_10-00-28_Cel_GAVsL4/tiny-deseq")
deseq_gaVl4<- read.csv("tinyrna_2023-07-21_10-00-28_cond1_Cel_WT_GA_cond2_Cel_WT_L4_deseq_table.csv")
df_up<-subset(deseq_gaVl4,deseq_gaVl4$log2FoldChange>2 & deseq_gaVl4$padj<0.05)
#2290 elements unique to GA
df_up<-subset(df_up,df_up$Classifier%in% c("Genes 26G ", "Pseudogene 26G", "Transposon 26G ", "miRNA 26G"))
df_down<-subset(deseq_gaVl4,deseq_gaVl4$log2FoldChange<(-2) & deseq_gaVl4$padj<0.05)
#4661 elements unique to L4

setwd("X:/Shamitha/tinyrna/START_HERE/tinyrna_2023-07-21_10-21-40_cbr_GAvL4/tiny-deseq")
deseq_gaVl4<- read.csv("tinyrna_2023-07-21_10-21-40_cond1_WT_L4_cond2_WT_GA_deseq_table.csv")

#L4 specific 26GRNA targets
df_up<-subset(deseq_gaVl4,deseq_gaVl4$log2FoldChange>2 & deseq_gaVl4$padj<0.05)
df_up<-subset(df_up,df_up$Classifier%in% c("Genes 26G ", "Pseudogene 26G", "Transposon 26G ", "miRNA 26G"))
write.csv(df_up, "X:/Shamitha/tinyrna/START_HERE/analysis/26g_targets_cbg_L4.csv")
#GA specific 26GRNA targets
df_down<-subset(deseq_gaVl4,deseq_gaVl4$log2FoldChange<(-2) & deseq_gaVl4$padj<0.05)
df_down<-subset(df_down,df_down$Classifier%in% c("Genes 26G ", "Pseudogene 26G", "Transposon 26G ", "miRNA 26G"))
write.csv(df_down, "X:/Shamitha/tinyrna/START_HERE/analysis/26g_targets_cbg_GA.csv")
