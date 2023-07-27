library('rtracklayer')
library('dplyr')
library("reshape2")
library("ggplot2")
library("seqinr")
library("tidyverse")
matisse_palette <- c("#FED766", "#F07167", "#353535", "#009975", "#004E64")
custom_palette <- c("#FFCBA4", "#F7A35C", "#FF7F00", "#0072B2", "#00538A")

setwd("X:/Shamitha/tinyrna/START_HERE/analysis/elegans")
targets_cel_l4<-read.csv("gtsf1_targets_cel_l4.csv")
unique(targets_cel_l4$Classifier)

#Pie chart pretttyy
l4_pie<-data.frame( Value= c(nrow(targets_cel_l4[targets_cel_l4$Classifier%in% c("Genes 26G "),]),
                              nrow(targets_cel_l4[targets_cel_l4$Classifier%in% c("Pseudogene 26G"),]),
                              nrow(targets_cel_l4[targets_cel_l4$Classifier%in% c("Transposon 26G "),]),
                              nrow(targets_cel_l4[targets_cel_l4$Classifier%in% c("miRNA 26G"),]),
                             nrow(targets_cel_l4[targets_cel_l4$Classifier%in% c("piRNA 26G"),])
                             ),
                    Category= c("Protein-coding gene", "Pseudogene", "Transposon", "miRNA", "piRNA")
                    )
l4_pie$Category<- factor(l4_pie$Category, levels = c("Protein-coding gene", "Pseudogene", "Transposon", "miRNA", "piRNA"))

#Pie chart pretttyy

ggplot(l4_pie, aes(x = "", y = Value, fill = Category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  labs(title = "GTSF-1 targets in C. elegans (L4)", fill = "Category")+
  scale_fill_manual(values = c("Protein-coding gene" = "#007EA7", 
                                "Pseudogene" =  "#D88198",
                                "Transposon" = "#FFB84E",
                                "miRNA" = "#55D6BE",
                               "piRNA" ="#F07167"))+
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))+
  geom_text(aes(label = Value), position = position_stack(vjust = 0.5))

#######################################################################################


##for briggsae

setwd("X:/Shamitha/tinyrna/START_HERE/analysis")
targets_cbr_GA<-read.csv("gtsf1_targets_cbg_GA.csv")
unique(targets_cbr_GA$Classifier)
GA_pie<-data.frame( Value= c(nrow(targets_cbr_GA[targets_cbr_GA$Classifier%in% c("Genes 26G "),]),
                             nrow(targets_cbr_GA[targets_cbr_GA$Classifier%in% c("Pseudogene 26G"),]),
                             nrow(targets_cbr_GA[targets_cbr_GA$Classifier%in% c("Transposon 26G "),]),
                             nrow(targets_cbr_GA[targets_cbr_GA$Classifier%in% c("piRNA 26G"),])),
                    Category= c("Protein-coding gene", "Pseudogene", "Transposon", "piRNA")
                    )


GA_pie$Category<- factor(GA_pie$Category, levels = c("Protein-coding gene", "Pseudogene", "Transposon", "piRNA"))
ggplot(GA_pie, aes(x = "", y = Value, fill = Category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  labs(title = "GTSF-1 targets in C.briggsae (GA)", fill = "Category")+
  scale_fill_manual(values = c("Protein-coding gene" = "#007EA7", 
                               "Pseudogene" =  "#D88198",
                               "Transposon" = "#FFB84E",
                               #"miRNA" = "#55D6BE",
                               "piRNA" ="#F07167"))+
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))+
  geom_text(aes(label = Value), position = position_stack(vjust = 0.5))

