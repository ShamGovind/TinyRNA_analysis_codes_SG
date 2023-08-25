#Volcano plots

##Project 16 Cbr L4
highlight_genes <- c("gtsf-1",
                     "rrf-3",
                     "eri-5",
                     "tofu-6",
                     "drh-3",
                    "tost-1",
                     "rde-12",
                     "pid-3",
                     "erh-2",
                     "tost-1",
                     "pid-1",
                     "ife-3")  # Genes to label
setwd("M:/Data/Projects/All IP-MS data/20230712_Shamitha_Cbriggsae_L4stage/analysis/20230710_results")
l4_ip<- read.csv("volcano_data_GTSF1.IgG.csv")
ggplot(l4_ip, aes(x = log2FC, y = X.log10.p.value.)) +
  geom_point(color = "grey", alpha = 0.6, size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#353535") + 
  geom_vline(xintercept= 1.5,linetype = "dashed", color = "#353535") + 
  geom_vline(xintercept= -1.5,linetype = "dashed", color = "#353535") +
  labs(title = "Volcano Plot Project 16 Cbr L4",
       x = "log2FC",
       y = "-log10(P-Value)")+
  theme_bw()+
      geom_point(data = l4_ip[l4_ip$plot_label %in% highlight_genes, ],  color = "grey2", size = 2 )+
  geom_point(data= subset (l4_ip, l4_ip$log2FC>1.5 & l4_ip$X.log10.p.value.>1.3) ,  color = "#66CC33", size = 3 )+
  geom_point(data= subset (l4_ip, l4_ip$log2FC< (-1.5) & l4_ip$X.log10.p.value.>1.3) ,  color = "#66CC33", size = 3 )+
geom_text(data = l4_ip[l4_ip$plot_label %in% highlight_genes, ], aes(label = plot_label), hjust = 0, vjust = 0)



##Plotting project 7
setwd( "U:/Data/Projects/All IP-MS data/Project 7/20200327_Project 7_20")
p7<- read.csv("volcano_proteinGroups.csv")
highlight_genes <- c("Cbr-gtsf-1",
                     "Cbr-rrf-3",
                     "Cbr-eri-5",
                     "Cbr-tofu-6",
                     "Cbr-drh-3",
                     "Cbr-tost-1",
                     "Cbr-rde-12",
                     "Cbr-pid-3",
                     "Cbr-erh-2",
                     "Cbr-tost-1",
                     "Cbr-pid-1",
                     "Cbr-ife-3")  # Genes to label
colnames(p7)
ggplot(p7, aes(x = difference_aCbrGTSF1_aIgG, y = log10.pvalue_aCbrGTSF1_aIgG)) +
  geom_point(color = "grey", alpha = 0.6, size = 2) +
  geom_hline(yintercept = -log10(0.05),linetype = "dashed", color = "#353535") + 
  geom_vline(xintercept= 1.5,linetype = "dashed", color = "#353535") + 
  geom_vline(xintercept= -1.5,linetype = "dashed", color = "#353535") + 
  labs(title = "Volcano Plot Project 6- Cbr WT",
       x = "Fold Change",
       y = "-log10(P-Value)")+
  xlim(-8,8)+
  theme_minimal()+
   geom_point(data = p7[p7$my_label %in% highlight_genes, ],  color = "grey2", size = 2 )+
  geom_point(data= subset (p7, p7$difference_aCbrGTSF1_aIgG>1.5 & p7$log10.pvalue_aCbrGTSF1_aIgG>1.30103) ,  color = "#66CC33", size = 3 )+
  geom_point(data= subset (p7, p7$difference_aCbrGTSF1_aIgG< (-1.5) & p7$log10.pvalue_aCbrGTSF1_aIgG>1.30103) ,  color = "#66CC33", size = 3 )+
 geom_text(data = p7[p7$my_label %in% highlight_genes, ], aes(label = my_label), hjust = 0, vjust = 0) 


ggplot(p7, aes(x = difference_aCbrGTSF1_aCelGTSF1, y = log10.pvalue_aCbrGTSF1_aCelGTSF1)) +
  geom_point(color = "grey", alpha = 0.6, size = 2) +
  geom_hline(yintercept = -log10(0.05),linetype = "dashed", color = "#353535") + # Add significance threshold (adjust if needed)
  geom_vline(xintercept= 1.5,linetype = "dashed", color = "#353535") + # Add significance threshold (adjust if needed)
  geom_vline(xintercept= -1.5,linetype = "dashed", color = "#353535") + 
  labs(title = "Volcano Plot Project 6- Cbr WT- CbrGTSF1 vs CelGTSF1",
       x = "Fold Change",
       y = "-log10(P-Value)")+
  xlim(-8,8)+
  theme_minimal()+
  geom_point(data = p7[p7$my_label %in% highlight_genes, ],  color = "grey2", size = 2 )+
  geom_point(data= subset (p7, p7$difference_aCbrGTSF1_aCelGTSF1>1.5 & p7$log10.pvalue_aCbrGTSF1_aCelGTSF1>1.30103) ,  color = "#66CC33", size = 3 )+
  geom_point(data= subset (p7, p7$difference_aCbrGTSF1_aCelGTSF1< (-1.5) & p7$log10.pvalue_aCbrGTSF1_aCelGTSF1>1.30103) ,  color = "#66CC33", size = 3 )+
  geom_text(data = p7[p7$my_label %in% highlight_genes, ], aes(label = my_label), hjust = 0, vjust = 0) 



##Project 14 ERI5 deletion
setwd( "U:/Data/Projects/All IP-MS data/20230712_Shamitha_Cbriggsae_ERI5strain/analysis/20230710_results")
p14<- read.csv("volcano_data_GTSF1.IgG.csv")
colnames(p14)
highlight_genes <- c("gtsf-1",
                     "rrf-3",
                     "eri-5",
                     "tofu-6",
                     "drh-3",
                     "tost-1",
                     "rde-12",
                     "pid-3",
                     "erh-2",
                     "tost-1",
                     "pid-1",
                     "ife-3")  # Genes to label

ggplot(p14, aes(x = log2FC, y = X.log10.p.value.)) +
  geom_point(color = "grey", alpha = 0.6, size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#353535") + 
  geom_vline(xintercept= 1.5,linetype = "dashed", color = "#353535") + 
  geom_vline(xintercept= -1.5,linetype = "dashed", color = "#353535") +
  labs(title = "Volcano Plot Project 14 Cbr eri5 del",
       x = "Fold Change",
       y = "-log10(P-Value)")+
  theme_bw()+
  xlim(-8,8)+
  geom_point(data = p14[p14$plot_label %in% highlight_genes, ],  color = "grey2", size = 2 )+
  geom_point(data= subset (p14, p14$log2FC>1.5 & p14$X.log10.p.value.>1.3) ,  color = "#66CC33", size = 3 )+
  geom_point(data= subset (p14, p14$log2FC< (-1.5) & p14$X.log10.p.value.>1.3) ,  color = "#66CC33", size = 3 )+
  geom_text(data = p14[p14$plot_label %in% highlight_genes, ], aes(label = plot_label), hjust = 0, vjust = 0)+
annotate("text", x=5, y =0, label ="C. briggsae ERI-5 ??; ??- GTSF-1", size = 6)+
  annotate("text", x= (-5), y =0, label ="C. briggsae ERI-5 ??; ??- IgG", size = 6)




#Project 15 Cbr GTSF1 del
setwd( "U:/Data/Projects/All IP-MS data/20230712_Shamitha_Cbriggsae_Control/analysis/20230710_results")
p15<- read.csv("volcano_data_GTSF1.IgG.csv")
colnames(p15)
highlight_genes <- c("gtsf-1",
                     "rrf-3",
                     "eri-5",
                     "tofu-6",
                     "drh-3",
                     "tost-1",
                     "rde-12",
                     "pid-3",
                     "erh-2",
                     "tost-1",
                     "pid-1",
                     "ife-3")  # Genes to label

ggplot(p15, aes(x = log2FC, y = X.log10.p.value.)) +
  geom_point(color = "grey", alpha = 0.6, size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#353535") + 
  geom_vline(xintercept= 1.5,linetype = "dashed", color = "#353535") + 
  geom_vline(xintercept= -1.5,linetype = "dashed", color = "#353535") +
  labs(title = "Volcano Plot Project 15- Cbr gtsf1 del",
       x = "Fold Change",
       y = "-log10(P-Value)")+
  theme_bw()+
  xlim(-8,8)+
  geom_point(data = p15[p15$plot_label %in% highlight_genes, ],  color = "grey2", size = 2 )+
  geom_point(data= subset (p15, p15$log2FC>1.5 & p15$X.log10.p.value.>1.3) ,  color = "#66CC33", size = 3 )+
  geom_point(data= subset (p15, p15$log2FC< (-1.5) & p15$X.log10.p.value.>1.3) ,  color = "#66CC33", size = 3 )+
  geom_text(data = p15[p15$plot_label %in% highlight_genes, ], aes(label = plot_label), hjust = 0, vjust = 0)+
  annotate("text", x=5, y =0, label ="C. briggsae WT", size = 6)+
  annotate("text", x= (-5), y =0, label ="C. briggsae GTSF-1 ?? ", size = 6)
  
  
  
data_15= subset (p15, p15$log2FC>1.5 & p15$X.log10.p.value.>1.3)
data_15<- data_15[, c(2,22,23)]




##Project12 PID3 IP
setwd("U:/Data/Projects/All IP-MS data/20230524_shamita_pid3_results/analysis/20230523_results")
p12<- read.csv("volcano_data_pid3.IgG.csv")
colnames(p12)
highlight_genes <- c("gtsf-1",
                     "rrf-3",
                     "eri-5",
                     "tofu-6",
                     "drh-3",
                     "tost-1",
                     "rde-12",
                     "pid-3",
                     "erh-2",
                     "tost-1",
                     "pid-1",
                     "ife-3")  # Genes to label

ggplot(p12, aes(x = log2FC, y = X.log10.p.value.)) +
  geom_point(color = "grey", alpha = 0.6, size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#353535") + 
  geom_vline(xintercept= 1.5,linetype = "dashed", color = "#353535") + 
  geom_vline(xintercept= -1.5,linetype = "dashed", color = "#353535") +
  labs(title = "Volcano Plot Project 12- Cbr PID3 IP",
       x = "Fold Change",
       y = "-log10(P-Value)")+
  theme_bw()+
  #xlim(-8,8)+
  geom_point(data = p12[p12$plot_label %in% highlight_genes, ],  color = "grey2", size = 2 )+
  geom_point(data= subset (p12, p12$log2FC>1.5 & p12$X.log10.p.value.>1.3) ,  color = "#66CC33", size = 3 )+
  geom_point(data= subset (p12, p12$log2FC< (-1.5) & p12$X.log10.p.value.>1.3) ,  color = "#66CC33", size = 3 )+
  geom_text(data = p12[p12$plot_label %in% highlight_genes, ], aes(label = plot_label), hjust = 0, vjust = 0)+
  annotate("text", x=5, y =0, label ="??- PID-3", size = 6)+
  annotate("text", x= (-5), y =0, label ="??- IgG", size = 6)
            

data_12= subset (p12, p12$log2FC>1.5 & p12$X.log10.p.value.>1.3)
data_12<- data_12[, c(2,28,29)]


data<- merge(data_12, data_15, by = "plot_label")


##Project13 Ppacificus IP 
library(ggplot2)
setwd("U:/Data/Projects/All IP-MS data/20230712_Shamitha_Ppacificus_GTSF1IP/analysis/20230710_results")
p13<- read.csv("volcano_data_WT.noGTSF1.csv")
colnames(p13)
highlight_genes <- c("gtsf-1",
                     "rrf-3",
                     "ekl-1",
                     "NYN domain protein",
                     "pir-1")  # Genes to label

ggplot(p13, aes(x = log2FC, y = X.log10.p.value.)) +
  geom_point(color = "grey", alpha = 0.6, size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#353535") + 
  geom_vline(xintercept= 1.5,linetype = "dashed", color = "#353535") + 
  geom_vline(xintercept= -1.5,linetype = "dashed", color = "#353535") +
  labs(title = "Project 13- GTSF-1 IP in P.pacificus Gravid-Adults",
       x = "Fold Change",
       y = "-log10(P-Value)")+
  theme_bw()+
  xlim(-8,8)+
  geom_point(data = p13[p13$plot_label %in% highlight_genes, ],  color = "grey2", size = 2 )+
  geom_point(data= subset (p13, p13$log2FC>1.5 & p13$X.log10.p.value.>1.3) ,  color = "#66CC33", size = 3 )+
  geom_point(data= subset (p13, p13$log2FC< (-1.5) & p13$X.log10.p.value.>1) ,  color = "#66CC33", size = 3 )+
  geom_text(data= subset (p13, p13$log2FC>1 & p13$X.log10.p.value.>1.3), aes(label = Protein.names),color="grey2" ,hjust = 0, vjust = 0)+
  geom_text(data = p13[p13$Protein.names %in% highlight_genes, ], aes(label = Protein.names), color= "blue", hjust = 0, vjust = 0)+
  annotate("text", x=5, y =0, label ="P. pacificus WT", size = 6)+
  annotate("text", x= (-5), y =0, label ="P. pacificus GTSF-1 del", size = 6)



##Plotting project 10- Ppacifius IP number 1
setwd( "U:/Data/Projects/All IP-MS data/20220528_Project10")
p10<- read.csv("volcano_proteinGroups.csv")
highlight_genes <- c("Cbr-gtsf-1",
                     "Cbr-rrf-3",
                     "Cbr-eri-5",
                     "Cbr-tofu-6",
                     "Cbr-drh-3",
                     "Cbr-tost-1",
                     "Cbr-rde-12",
                     "Cbr-pid-3",
                     "Cbr-erh-2",
                     "Cbr-tost-1",
                     "Cbr-pid-1",
                     "Cbr-ife-3")  # Genes to label
colnames(p7)
ggplot(p7, aes(x = difference_aCbrGTSF1_aIgG, y = log10.pvalue_aCbrGTSF1_aIgG)) +
  geom_point(color = "grey", alpha = 0.6, size = 2) +
  geom_hline(yintercept = -log10(0.05),linetype = "dashed", color = "#353535") + 
  geom_vline(xintercept= 1.5,linetype = "dashed", color = "#353535") + 
  geom_vline(xintercept= -1.5,linetype = "dashed", color = "#353535") + 
  labs(title = "Volcano Plot Project 6- Cbr WT",
       x = "Fold Change",
       y = "-log10(P-Value)")+
  xlim(-8,8)+
  theme_minimal()+
  geom_point(data = p7[p7$my_label %in% highlight_genes, ],  color = "grey2", size = 2 )+
  geom_point(data= subset (p7, p7$difference_aCbrGTSF1_aIgG>1.5 & p7$log10.pvalue_aCbrGTSF1_aIgG>1.30103) ,  color = "#66CC33", size = 3 )+
  geom_point(data= subset (p7, p7$difference_aCbrGTSF1_aIgG< (-1.5) & p7$log10.pvalue_aCbrGTSF1_aIgG>1.30103) ,  color = "#66CC33", size = 3 )+
  geom_text(data = p7[p7$my_label %in% highlight_genes, ], aes(label = my_label), hjust = 0, vjust = 0) 


ggplot(p7, aes(x = difference_aCbrGTSF1_aCelGTSF1, y = log10.pvalue_aCbrGTSF1_aCelGTSF1)) +
  geom_point(color = "grey", alpha = 0.6, size = 2) +
  geom_hline(yintercept = -log10(0.05),linetype = "dashed", color = "#353535") + # Add significance threshold (adjust if needed)
  geom_vline(xintercept= 1.5,linetype = "dashed", color = "#353535") + # Add significance threshold (adjust if needed)
  geom_vline(xintercept= -1.5,linetype = "dashed", color = "#353535") + 
  labs(title = "Volcano Plot Project 6- Cbr WT- CbrGTSF1 vs CelGTSF1",
       x = "Fold Change",
       y = "-log10(P-Value)")+
  xlim(-8,8)+
  theme_minimal()+
  geom_point(data = p7[p7$my_label %in% highlight_genes, ],  color = "grey2", size = 2 )+
  geom_point(data= subset (p7, p7$difference_aCbrGTSF1_aCelGTSF1>1.5 & p7$log10.pvalue_aCbrGTSF1_aCelGTSF1>1.30103) ,  color = "#66CC33", size = 3 )+
  geom_point(data= subset (p7, p7$difference_aCbrGTSF1_aCelGTSF1< (-1.5) & p7$log10.pvalue_aCbrGTSF1_aCelGTSF1>1.30103) ,  color = "#66CC33", size = 3 )+
  geom_text(data = p7[p7$my_label %in% highlight_genes, ], aes(label = my_label), hjust = 0, vjust = 0) 



unique(colnames(p10))
p10_hits<- p10_hits [, c(98,117,133,134)]
p10_hits<- subset(p10, p10$log10.pvalue_aPpGTSF1_aIgG>1 & p10$difference_aPpGTSF1_aIgG>1)

unique(colnames(p13))
p13_hits<- p13_hits[ , c(2, 28,29,35)]
p13_hits<- subset(p13, p13$log2FC>1 & p13$X.log10.p.value.>1)

ppa_hits<- subset(p10_hits, p10_hits$WB_ID %in% p13_hits$plot_label)
ppa_hits<- subset(p13_hits,  p13_hits$plot_label %in% p10_hits$WB_ID)
