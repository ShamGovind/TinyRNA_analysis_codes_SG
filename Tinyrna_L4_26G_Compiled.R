#06-08-23 Let's do the following for L4 samples:
#Define targets with WTvsDel plot 
#Define SetABC for both- Boxplot, MA plot, Scatterplot between species
#Scatterplot to compare SetC

#For L4 and GA- Plot 26G/gene distribution


library(ggplot2)
library(dplyr)
setwd("X:/imb-kettinggr/Shamitha/tinyrna/START_HERE/tinyrna_2023-07-17_21-29-34_Cel_2021_08_L4_26G/tiny-deseq")
ncount_cel_og <- read.csv("tinyrna_2023-07-17_21-29-34_norm_counts.csv")
unique(ncount_cel_og$Classifier)
ncount_cel<- subset(ncount_cel_og, ncount_cel_og$Classifier%in% c("Genes 26G ", "Pseudogene 26G", "Transposon 26G", "piRNA 26G", "miRNA 26G"))
ncount_cel<- data.frame(Gene_elegans = ncount_cel$Feature.ID,
                        Classifier_elegans = ncount_cel$Classifier,
                        Feature_elegans= ncount_cel$Feature.Name,
                        Count_del= (ncount_cel$Cel_Del_L4_rep_1 + ncount_cel$Cel_Del_L4_rep_2 + ncount_cel$Cel_Del_L4_rep_3)/3,
                        Count_wt= (ncount_cel$Cel_WT_L4_rep_1 + ncount_cel$Cel_WT_L4_rep_2 + ncount_cel$Cel_WT_L4_rep_3)/3)

setwd("X:/imb-kettinggr/Shamitha/tinyrna/START_HERE/tinyrna_2023-07-18_14-02-56_cbr_26G_L4/tiny-deseq")
ncount_cbr_og<- read.csv("tinyrna_2023-07-18_14-02-56_norm_counts.csv")
ncount_cbr<- subset(ncount_cbr_og, ncount_cbr_og$Classifier%in% c("Genes 26G ", "Pseudogene 26G", "Transposon 26G", "piRNA 26G", "miRNA 26G"))
ncount_cbr<- data.frame(Gene_briggsae = ncount_cbr$Feature.ID,
                        Classifier_briggsae = ncount_cbr$Classifier,
                        Feature_briggsae= ncount_cbr$Feature.Name,
                        Count_del= (ncount_cbr$gtsf.1.xf345._L4_rep_1 + ncount_cbr$gtsf.1.xf345._L4_rep_2 + ncount_cbr$gtsf.1.xf345._L4_rep_3)/3,
                        Count_wt= (ncount_cbr$WT_L4_rep_1 + ncount_cbr$WT_L4_rep_2 + ncount_cbr$WT_L4_rep_3)/3)

#use te next two commands for plotting only 
ncount_cbr<- subset(ncount_cbr, ncount_cbr$Count_del + ncount_cbr$Count_wt> 0)
ncount_cel<- subset(ncount_cel, ncount_cel$Count_del + ncount_cel$Count_wt> 0)
ggplot(ncount_cel, aes(y = log2(Count_del), x = log2(Count_wt))) +
  geom_point(color = "#353535", size = 1) +
  geom_abline(color= "grey")+
  theme_minimal()+
  #xlim(c(0,10))+
  #ylim(c(0,10))+
  labs(y = "Log2 (Mean normalized counts GTSF1 del)",  x = " Log2 (Mean normalized counts WT)", title = "26G counts per gene in elegans L4 animals", subtitle =  "Targets= 2417/6306")+
  geom_point(data = subset(ncount_cel, ncount_cel$Count_wt/ncount_cel$Count_del>2), color = "#F07167", size = 1)+
  geom_point(data = subset(ncount_cel, ncount_cel$Count_wt<5 & ncount_cel$Count_del<5), color = "grey", size = 1)


nrow(subset(ncount_cel, ncount_cel$Count_wt>5 &  ncount_cel$Count_wt/ncount_cel$Count_del>2))

#Bar plot of target biotypes
targets_cel<- subset(ncount_cel, ncount_cel$Count_wt>1 &  ncount_cel$Count_wt/ncount_cel$Count_del>2)
targets_cbr<- subset(ncount_cbr, ncount_cbr$Count_wt>1 &  ncount_cbr$Count_wt/ncount_cbr$Count_del>2)

unique(targets_cel$Classifier_elegans)
unique(targets_cbr$Classifier_briggsae)

biotype<- data.frame(Psuedogene_cel<- nrow(subset(targets_cel, targets_cel$Classifier_elegans=="Pseudogene 26G"))/ nrow(targets_cel),
                     Psuedogene_cbr<- nrow(subset(targets_cbr, targets_cbr$Classifier_briggsae=="Pseudogene 26G")) / nrow(targets_cbr),
                     Gene_cel<- nrow(subset(targets_cel, targets_cel$Classifier_elegans=="Genes 26G "))/ nrow(targets_cel),
                     Gene_cbr<- nrow(subset(targets_cbr, targets_cbr$Classifier_briggsae=="Genes 26G ")) / nrow(targets_cbr))
library(reshape2)
biotype<- melt(biotype)
biotype$species<- c( "Cel", "Cbr", "Cel", "Cbr")
biotype$type<- rep(c("Psuedogene", "Protien-Coding gene"), each=2)

library(RColorBrewer)
ggplot(biotype, aes(x=species, y=value, fill =type))+
  geom_bar(stat = "identity", position = "dodge", width = 0.4)+
  scale_fill_brewer(palette = "Dark2")+theme_bw()

nrow(subset(targets_cel, targets_cel$Classifier_elegans=="Pseudogene 26G"))
nrow(subset(targets_cbr, targets_cbr$Classifier_briggsae=="Pseudogene 26G"))
nrow(subset(targets_cel, targets_cel$Classifier_elegans=="Genes 26G "))
nrow(subset(targets_cbr, targets_cbr$Classifier_briggsae=="Genes 26G "))


##Scatterplots of counts/target between species and identification of homologs
setwd("X:/imb-kettinggr/Shamitha/BLAST/ncbi-blast-2.14.0+/blastdb")
homolog_master<- read.csv("CelVCbr_homologs_prot_gene.csv")
homolog_gene<- read.csv("CelVCbr_homologs_unique_gene_only.csv")
homolog_master<- homolog_master[,-1]
homolog_gene<- homolog_gene[,-1]

colnames(homolog_gene)<- c("Gene_elegans", "Gene_briggsae")
a<-merge(homolog_gene,ncount_cel, by= "Gene_elegans")
a<-merge(a,ncount_cbr, by = "Gene_briggsae")
colnames(a)<- c( "Gene_briggsae", "Gene_elegans", "Classifier_elegans", "Feature_elegans","Count_del.elegans","Count_wt.elegans","Classifier_briggsae", "Feature_briggsae","Count_del.briggsae","Count_wt.briggsae")
a<- subset(a, a$Count_wt.elegans + a$Count_wt.briggsae> 0)




#Categorize genes from matrix a into set B C, 1 IS for Celegans and 2 is for Cbriggsae:
ncount_cel<- subset(ncount_cel, ncount_cel$Count_del + ncount_cel$Count_wt> 0)
ncount_cel<- subset(ncount_cel,  ncount_cel$Count_wt/ncount_cel$Count_del > 2)
SetA1<- ncount_cel[!(ncount_cel$Gene_elegans%in%homolog_gene$Gene_elegans),]
SetA1<- subset(SetA1, SetA1$Count_wt>5)

#B1= Homologs+ targets in elegans(wt>1, wt/del>2) and not a target in briggsae (wt<1)
SetB1<- subset(a, a$Count_wt.elegans>0 & a$Count_wt.briggsae<5)
SetB1<- subset(SetB1, SetB1$Count_wt.elegans/SetB1$Count_del.elegans>2)  
SetB1<- subset(SetB1, SetB1$Count_wt.elegans>5)

ncount_cbr<- subset(ncount_cbr, ncount_cbr$Count_del + ncount_cbr$Count_wt> 0)
ncount_cbr<- subset(ncount_cbr,  ncount_cbr$Count_wt/ncount_cbr$Count_del > 2)
SetA2<- ncount_cbr[!(ncount_cbr$Gene_briggsae%in%homolog_gene$Gene_briggsae),]
SetA2<- subset(SetA2, SetA2$Count_wt>5)

SetB2<- subset(a, a$Count_wt.briggsae>0 & a$Count_wt.elegans<5)
SetB2<- subset(SetB2, SetB2$Count_wt.briggsae/SetB2$Count_del.briggsae>2) 
SetB2<- subset(SetB2, SetB2$Count_wt.briggsae>5)

SetC<- subset(a, (a$Count_wt.briggsae>5 & a$Count_wt.elegans>5))
SetC<- subset(SetC, SetC$Count_wt.briggsae/SetC$Count_del.briggsae>2 & SetC$Count_wt.elegans/SetC$Count_del.elegans>2)


x<-subset(SetB1, SetB1$Gene_elegans%in%SetC$Gene_elegans)

#Fisher's test
data1 <- matrix(c(nrow(homolog_gene), nrow(SetB1), nrow(SetB2), nrow(SetC)), nrow = 2, byrow = TRUE)
colnames(data1) <- c(" Homolog in elegans", "Target in elegans")
rownames(data1) <- c("Homolog in briggsae", "Target in briggsae")

result1 <- fisher.test(data1)
print(result)

#Fisher's test
data2 <- matrix(c(1718, 1294, 473, 2532), nrow = 2, byrow = TRUE)
colnames(data2) <- c( "Target in elegans","Target in briggsae")
rownames(data2) <- c( "Homolog in elegans","Homolog in briggsae")

result2 <- fisher.test(data2)
print(result2)




#Scatterplot of meancounts/genes between species for SetC genes
ggplot(SetC, aes(y = log2(Count_wt.elegans), x = log2(Count_wt.briggsae))) +
  geom_point(color = "#353535", size = 1) +
  geom_smooth(method = "lm", color = "grey")+
  theme_minimal()+
  xlim(c(2,10))+
  ylim(c(2,10))+
  geom_text(x = 8, y = 10, label ="n= 593, rho= 0.31, a= 2.9e-14", color = "#0072B2")+
  labs(y = "Log2 (GTSF1 WT counts elegans )",  x = " Log2 (GTSF1 WT counts briggsae)", title = "26G counts/gene : L4 SET C")

  cor.test(SetC$Count_wt.elegans, SetC$Count_wt.briggsae, method = "spearman")



  b<- subset(a, a$Count_wt.briggsae/a$Count_del.briggsae>2 & a$Count_wt.elegans/a$Count_del.elegans>2)
  b<- subset(b, !(b$Count_wt.elegans<5 & b$Count_wt.briggsae<5)) 
  ggplot(b, aes(y = log2(Count_wt.elegans), x = log2(Count_wt.briggsae))) +
    geom_point(color = "grey", size = 1.5) +
    geom_abline(color= "grey")+
    theme_minimal()+
    geom_point(data = subset(b, b$Gene_elegans %in% SetB1$Gene_elegans), color = "royalblue", size = 1.5)+
    geom_point(data = subset(b, b$Gene_briggsae %in% SetB2$Gene_briggsae), color = "#0072B2", size = 1.5)+
    geom_point(data = subset(b, b$Gene_briggsae %in% SetC$Gene_briggsae), color = "#F07167", size = 1.5)+
    geom_point(data = subset(b, b$Gene_elegans %in% SetC$Gene_elegans), color = "#F07167", size = 1.5)+
     labs(y = "Log2 (WT counts elegans )",  x = " Log2 (WT counts briggsae)", title = "26G counts/gene : L4")
  #geom_point(data = subset(a, a$Count_del.elegans+a$Count_del.briggsae<1), color = "#F07167", size = 1)
  

  #GO profiling of these sets
  library(clusterProfiler)
  library("org.Ce.eg.db")
  
  #Create the background for GO analysis:
  #For A1: All annotated genes
  #For B1, C: All homologous genes
  
  go_result_B2 <- enrichGO(
    SetB2$Gene_elegans,
    OrgDb = org.Ce.eg.db,
    keyType = "WORMBASE",
    #universe=ncount_cel$Gene_elegans,
    ont = "BP",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    readable = TRUE,
    pool = FALSE
  )
  
enrichplot::dotplot(go_result_B2)

library("ggplot2")
setwd("X:/imb-kettinggr/Shamitha/tinyrna/START_HERE/tinyrna_2023-07-18_14-02-56_cbr_26G_L4/tiny-deseq")
deseq_l4<- read.csv("tinyrna_2023-07-18_14-02-56_cond1_gtsf-1(xf345)_L4_cond2_WT_L4_deseq_table.CSV")
ggplot(deseq_l4, aes(y = log2FoldChange, x = log2(baseMean))) +
  geom_point(color = "grey", size = 0.5) +
  theme_minimal() +
  labs(y = "Log2 Fold Change", x = "Log2 Mean of Counts", title = "26G RNA targets in C. briggsae L4") +
  #geom_point(data = subset(deseq_l4, deseq_l4$padj<0.01), color = "#F07167", size = 1)+
  geom_point(data = subset(deseq_l4, deseq_l4$Feature.ID %in% SetA2$Gene_briggsae), color = "#009975", size = 1)+
  geom_point(data = subset(deseq_l4, deseq_l4$Feature.ID %in% SetC$Gene_briggsae), color = "#F07167", size = 1)+
  geom_point(data = subset(deseq_l4, deseq_l4$Feature.ID %in% SetB2$Gene_briggsae), color = "#0072B2", size = 1)

df<- data.frame(SetA2= length(SetA2$Gene_briggsae),
                SetB2= length(SetB2$Gene_briggsae),
                SetC= length(SetC$Gene_briggsae))

library(reshape2)
df<- melt(df)

df$variable <- as.factor(df$variable)
levels(df$variable) <- c("SetA1", "SetA2", "SetB1", "SetB2", "SetC")

df$group<- c("Novel targets", "Homologs + unique targets" , "Novel targets", "Homologs + unique targets" ,"Homologs + Shared targets")
ggplot(df, aes(x = variable, y = value, fill=variable)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  ylim(c(0,1300))+
  scale_fill_manual(values = c("#009975","#0072B2","#F07167"))+
  #labs(title = "Grouped Bar Chart",x = "Category",y = "Value",fill = "Group") +
  theme_minimal()
