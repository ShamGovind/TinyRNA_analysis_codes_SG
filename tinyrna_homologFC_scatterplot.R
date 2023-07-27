#Compare 26G RNA targets with homologs table
library(rtracklayer)
library(seqinr)
library(ggplot2)



setwd("X:/Shamitha/BLAST/ncbi-blast-2.14.0+/blastdb")
homolog_master<- read.csv("CelVCbr_homologs_prot_gene.csv")
homolog_gene<- read.csv("CelVCbr_homologs_unique_gene_only.csv")
homolog_master<- homolog_master[,-1]
homolog_gene<- homolog_gene[,-1]

setwd("X:/Shamitha/tinyrna/START_HERE/analysis/elegans")
targets_cel_l4<-read.csv("gtsf1_targets_cel_L4_hc.csv")
targets_cel_ga<-read.csv("gtsf1_targets_cel_ga_hc.csv")

setwd("X:/Shamitha/tinyrna/START_HERE/analysis/briggsae")
targets_cbg_l4<-read.csv("gtsf1_targets_cbg_L4_hc.csv")
targets_cbg_ga<-read.csv("gtsf1_targets_cbg_GA_hc.csv")



#######
setwd("X:/Shamitha/tinyrna/START_HERE/results/tinyrna_2023-07-17_20-17-19_Cel_2021_08_GA_26G/tiny-deseq")
deseq_ga_cel<- read.csv("tinyrna_2023-07-17_20-17-19_cond1_Cel_Del_GA_cond2_Cel_WT_GA_deseq_table.csv")
setwd("X:/Shamitha/tinyrna/START_HERE/tinyrna_2023-07-18_15-19-45_cbr_26G_GA/tiny-deseq")
deseq_ga_cbr<- read.csv("tinyrna_2023-07-18_15-19-45_cond1_gtsf-1(xf345)_GA_cond2_WT_GA_deseq_table.csv")


#Plot log2fc of GA homologs
hom_ga_cel<- deseq_ga_cel[deseq_ga_cel$Feature.ID%in%homolog_gene$Gene_Celegans,c(1,2,5)]
hom_ga_cbr<- deseq_ga_cbr[deseq_ga_cbr$Feature.ID%in%homolog_gene$Gene_Cbrigsae,c(1,2,5)]
colnames(hom_ga_cel)<- c("Gene_Celegans", "Classifier_elegans", "log2FC_elegans")
colnames(hom_ga_cbr)<- c("Gene_Cbrigsae", "Classifier_briggsae", "log2FC_briggsae")
a<-merge(homolog_gene,hom_ga_cel, by= "Gene_Celegans")
a<-merge(a,hom_ga_cbr, by = "Gene_Cbrigsae")
?merge
a <- replace(a, is.na(a), 0)
#list of common targets (conserved targets)
b<-subset(a, a$Gene_Celegans %in% targets_cel_ga$Feature.ID)
b<-subset(b, b$Gene_Cbrigsae%in%targets_cbg_ga$Feature.ID)


ggplot(a, aes(y = log2FC_elegans, x = log2FC_briggsae)) +
  geom_point(color = "grey", size = 1) +
  geom_abline(color= "#353535")+
  theme_minimal() +
  labs(y = "log2FC elegans",  x = "log2FC briggsae", title = "26G fold-change of all homologs genes in GA animals", subtitle= "26G-RNA target genes highlighted
       log2FC=0 indicates NA values")+
  xlim(c(-11,5))+
  ylim(c(-11,5))+
  geom_point(data = subset(a, a$Gene_Celegans %in% targets_cel_ga$Feature.ID), color = "#009975", size = 1.2)+
  geom_point(data = subset(a, a$Gene_Cbrigsae %in% targets_cbg_ga$Feature.ID), color = "#0072B2", size = 1.2)+
  geom_point(data=b, color="#F07167",  size = 1.2)
  

setwd("X:/Shamitha/tinyrna/START_HERE/results/tinyrna_2023-07-17_21-29-34_Cel_2021_08_L4_26G/tiny-deseq")
deseq_l4_cel<- read.csv("tinyrna_2023-07-17_21-29-34_cond1_Cel_Del_L4_cond2_Cel_WT_L4_deseq_table.csv")
setwd("X:/Shamitha/tinyrna/START_HERE/tinyrna_2023-07-18_14-02-56_cbr_26G_L4/tiny-deseq")
deseq_l4_cbr<- read.csv("tinyrna_2023-07-18_14-02-56_cond1_gtsf-1(xf345)_L4_cond2_WT_L4_deseq_table.csv")

#Plot log2FC of L4 homologs

hom_l4_cel<- deseq_l4_cel[deseq_l4_cel$Feature.ID%in%homolog_gene$Gene_Celegans,c(1,2,5)]
hom_l4_cbr<- deseq_l4_cbr[deseq_l4_cbr$Feature.ID%in%homolog_gene$Gene_Cbrigsae,c(1,2,5)]
colnames(hom_l4_cel)<- c("Gene_Celegans", "Classifier_elegans", "log2FC_elegans")
colnames(hom_l4_cbr)<- c("Gene_Cbrigsae", "Classifier_briggsae", "log2FC_briggsae")
a<-merge(homolog_gene,hom_l4_cel, by= "Gene_Celegans")
a<-merge(a,hom_l4_cbr, by = "Gene_Cbrigsae")
?merge
a <- replace(a, is.na(a), 0)
#list of common targets (conserved targets)
b<-subset(a, a$Gene_Celegans %in% targets_cel_l4$Feature.ID)
b<-subset(b, b$Gene_Cbrigsae%in%targets_cbg_l4$Feature.ID)


ggplot(a, aes(y = log2FC_elegans, x = log2FC_briggsae)) +
  geom_point(color = "grey", size = 00.8) +
  geom_abline(color= "#353535")+
  theme_minimal() +
  labs(y = "log2FC elegans",  x = "log2FC briggsae", title = "26G fold-change of all homologs genes in L4 animals", subtitle= "26G-RNA target genes highlighted
       log2FC=0 indicates NA values")+
  xlim(c(-12,5))+
  ylim(c(-12,5))+
  geom_point(data = subset(a, a$Gene_Celegans %in% targets_cel_l4$Feature.ID), color = "#009975", size = 0.8)+
  geom_point(data = subset(a, a$Gene_Cbrigsae %in% targets_cbg_l4$Feature.ID), color = "#0072B2", size = 0.8)+
  geom_point(data=b, color="#F07167",  size = 0.8)


###########
#plot 26g rna levels for

BaseMean_cel_ga=targets_cel_ga$baseMean
BaseMean_cel_l4=targets_cel_l4$baseMean
BaseMean_cbr_ga=targets_cbg_ga$baseMean
BaseMean_cbr_l4=targets_cbg_l4$baseMean
               
data <- data.frame(
  Value = c(BaseMean_cel_ga, BaseMean_cbr_ga, BaseMean_cel_l4, BaseMean_cbr_l4),
  Group = factor(c (rep("GA elegans ", length(BaseMean_cel_ga)), rep("GA briggsae", length(BaseMean_cbr_ga)), rep("L4 elegans", length(BaseMean_cel_l4)), rep("L4 briggsae", length(BaseMean_cbr_l4)) )),
  Stage = c(rep("GA", length(BaseMean_cel_ga)+ length(BaseMean_cbr_ga)), rep("L4", length(BaseMean_cel_l4)+ length(BaseMean_cbr_l4)))
  )

# Create the boxplot with data points

ggplot(data, aes(x = Group, y = Value)) +
  geom_boxplot(outlier.shape = NA, color= "grey") +
  geom_point(data = subset(data, data$Stage== "GA"), position = position_jitter(width = 0.2), color = "#009975", size = 0.8)+
  geom_point(data = subset(data, data$Stage== "L4"), position = position_jitter(width = 0.2), color = "#0072B2", size = 0.8)+
  #geom_point(data = subset(data, data$Value>1000),  fill = "red", shape = 24, size = 3) +
  labs(x= NULL, y = "26G RNA BaseMean") +
  ylim(c(-100,2000))+
  ggtitle("26G RNA BaseMean per target gene") +
  theme_minimal()


