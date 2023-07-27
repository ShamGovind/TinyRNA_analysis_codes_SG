#Compare 26G RNA targets with homologs table
library(rtracklayer)
library(seqinr)
library(ggplot2)
install.packages("VennDiagram")
library(VennDiagram)

#Load appropriate files

setwd("X:/Shamitha/BLAST/ncbi-blast-2.14.0+/blastdb")
homolog_master<- read.csv("CelVCbr_homologs_prot_gene.csv")
homolog_gene<- read.csv("CelVCbr_homologs_unique_gene_only.csv")
homolog_master<- homolog_master[,-1]
homolog_gene<- homolog_gene[,-1]


#Prior to loading, I modified the targets csv files on excel to remove the Gene: character- also removed first column which was rownames. 
setwd("X:/Shamitha/tinyrna/START_HERE/analysis")
targets_cel_l4<-read.csv("gtsf1_targets_cel_L4_hc.csv")
targets_cel_ga<-read.csv("gtsf1_targets_cel_ga_hc.csv")

setwd("X:/Shamitha/tinyrna/START_HERE/analysis")
targets_cbg_l4<-read.csv("gtsf1_targets_cbg_L4_hc.csv")
targets_cbg_ga<-read.csv("gtsf1_targets_cbg_GA_hc.csv")


#Q1:L4 stage has 4-10 times more target genes than GA stage. Are there any unique targets in GA?
nrow(targets_cel_ga[!(targets_cel_ga$Feature.ID%in%targets_cel_l4$Feature.ID),])
nrow(targets_cbg_ga[!(targets_cbg_ga$Feature.ID%in%targets_cbg_l4$Feature.ID),])

#A:Hmm, In celegans, there are 16 26G RNA targets unique to GA AND in Cbriggsae 90!!! 


#VennDiagram:
GA = targets_cbg_ga$Feature.ID
L4 = targets_cbg_l4$Feature.ID
# Load the required library
library(VennDiagram)

# Create the Venn diagram
venn.diagram(
  x = list(Set1 = GA, Set2 = L4),
  filename = "venn_cbrigsae.png",  # Set filename to NULL for interactive display
  category.names = c("GA", "L4"),
  output = TRUE,
  imagetype = "png",
  fill = c("#009975", "#0072B2"),
  alpha = 0.5,
  label.col = "black",
  cex = 2.5,
  fontfamily = "sans",
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = 0,
  margin = 0.2,
  main = "Overlap of 26G RNA targets in C. briggsae",
  main.cex = 2,
  title = "Sample Title",
  title.col = "black",
  title.cex = 1.8
)



library(dplyr)
#Q2a: Which 26G target genes from C. briggsae have a homolog in C. elegans?
colnames(targets_cbg_ga)<- c("Gene_Cbrigsae","Classifier", "Feature.Name","baseMean","log2FoldChange", "lfcSE","stat","pvalue","padj")
cbgGA_hom<- merge(homolog_gene,targets_cbg_ga, by="Gene_Cbrigsae")

colnames(targets_cbg_l4)<-c("Gene_Cbrigsae","Classifier", "Feature.Name","baseMean","log2FoldChange", "lfcSE","stat","pvalue","padj")
cbgL4_hom<- merge(homolog_gene,targets_cbg_l4, by="Gene_Cbrigsae")
#A: The conserved targets are saved in cbgGA and cbgL4 _hom tables

#Q2b What proportion of Cbr 26G protein coding target genes are conserved in Cel
length(unique(cbgGA_hom$Gene_Cbrigsae))/nrow(subset(targets_cbg_ga,targets_cbg_ga$Classifier =="Genes 26G "))
length(unique(cbgL4_hom$Gene_Cbrigsae))/nrow(subset(targets_cbg_l4,targets_cbg_l4$Classifier =="Genes 26G "))
#A:0.231 and 0.765
df<- data.frame(sample= c("Gravid Adults","L4s"), 
                  value= c(0.231, 0.765))

ggplot(df, aes(x=sample, y=value))+
  geom_bar(stat = "identity", width=0.7,
           fill = c("#009975", "#0072B2"))+
  labs(title = "26G RNA targets from L4 stage are more conserved than from GA stage",
       x= NULL,
       y= "Proportion of genes with a homolog in C. elegans")+ylim(0,1)+theme_minimal()
  #geom_text(aes(label = Value))+
  
                
                
#Q2c:Are the 75 conserved targets from Cbr GA samples present in L4 samples?
a<-targets_cbg_ga
b<-cbgGA_hom[(cbgGA_hom$Gene_Cbrigsae%in%targets_cbg_l4$Gene_Cbrigsae),]
c<-cbgGA_hom[!(cbgGA_hom$Gene_Cbrigsae%in%targets_cbg_l4$Gene_Cbrigsae),]
remove(a,b,c)

df<- data.frame(sample= c("All","Shared with L4", "Unique to GA"), 
                value= c(0.231, 0.231, 0))

ggplot(df, aes(x=sample, y=value))+
  geom_bar(stat = "identity", width=0.8,
           fill = c("#009975", "#32AA8C", "#005C48"))+
  labs(title = "All conserved target genes in GA samples are also targetted in L4 samples",
       x= "26G RNA target genes in GA",
       y= "Proportion of conserved genes")+ylim(0,1)+theme_minimal()
remove(df)

a$Gene_Cbrigsae%in%homolog_gene$Gene_Cbrigsae
#A: Yes, all of them. Thus, the targets unique to GA samples are not conserved in Celegans



#Q3: Is the homologous gene a 26G RNA target in C. elegans?

cbgGA_hom$Gene_Celegans%in%targets_cel_l4$Feature.ID
subset(cbgGA_hom,cbgGA_hom$Gene_Celegans%in%targets_cel_l4$Feature.ID)

a<-subset(cbgL4_hom,cbgL4_hom$Gene_Celegans%in%targets_cel_l4$Feature.ID)
a<-length(unique(a$Gene_Cbrigsae))
a2<-subset(cbgL4_hom,cbgL4_hom$Gene_Celegans%in%targets_cel_ga$Feature.ID)
a2<-length(unique(a2$Gene_Cbrigsae))

b<-subset(cbgGA_hom,cbgGA_hom$Gene_Celegans%in%targets_cel_ga$Feature.ID)
b<-length(unique(b$Gene_Cbrigsae))
b2<-subset(cbgGA_hom,cbgGA_hom$Gene_Celegans%in%targets_cel_l4$Feature.ID)
b2<-length(unique(b2$Gene_Cbrigsae))

#Answer, yes 527 out of 1794 L4 targets are also targeted by 26G RNAs in C.elegans 

x<-(nrow(targets_cbg_l4))-(length(unique(cbgL4_hom$Gene_Cbrigsae)))
y<-length(unique(cbgL4_hom$Gene_Cbrigsae))-a-a2
z<-a
l<-a2
q<-nrow(targets_cbg_l4)
  
x2<-(nrow(targets_cbg_ga))-(length(unique(cbgGA_hom$Gene_Cbrigsae)))
y2<-length(unique(cbgGA_hom$Gene_Cbrigsae))-b-b2
z2<-b
l2<-b2
q2<-(nrow(targets_cbg_ga))

df<- data.frame(label= c("No homolog", "homolog in C.elegans", "homolog+26G RNA target in C. elegans L4", "homolog+26G RNA target in C. elegans GA" ),
                L4=c(x/q,y/q,z/q,l/q),
                GA=c(x2/q2, y2/q2,l2/q2, z2/q2)
                )

library(reshape2)
library(RColorBrewer)
df<-melt(df)
ggplot(df, aes(x = variable, y = value, fill = label)) +
  scale_fill_manual(values = c("#F07167","#009975","#0072B2","#FFCBA4"))+
  geom_bar(stat = "identity", width = 0.6) +
  labs(title = "Conservation of C. briggsae 26G RNA target genes ",
       x = NULL,
       y = NULL,
       fill = "Target gene with:")+
  theme_minimal()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 14, face = "bold"))

targets_cel_ga[targets_cel_ga$Gene_Ce%in%cbgGA_hom$Gene_Celegans,]
nrow(targets_cel_l4[targets_cel_l4$Feature.ID%in%cbgGA_hom$Gene_Celegans,])
#Answer, yes 34 of these are also targeted by 26G RNAs in C.elegans - a set of conserved targets. 



#Wait, what if I look at the targets unique to GA?
nrow(targets_cbg_ga[(targets_cbg_ga$Gene_Cbrigsae%in%targets_cbg_l4$Gene_Cbrigsae),])/ nrow(targets_cbg_ga)
nrow(targets_cbg_ga[!(targets_cbg_ga$Gene_Cbrigsae%in%targets_cbg_l4$Gene_Cbrigsae),])
#None of these 71 genes have a cel homolog, they are unique to cbrigsae