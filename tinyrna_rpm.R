library('rtracklayer')
library('dplyr')
library("reshape2")
library("ggplot2")

#For Cbregans GA samples 26G and 21U reads only
setwd("X:/Shamitha/tinyrna/START_HERE/tinyrna_2023-07-17_20-17-19_Cel_2021_08_GA_26G/tiny-count")
Cel_GA_fc<-read.csv('tinyrna_2023-07-17_20-17-19_feature_counts.csv')
Cel_GA_stat<-read.csv('tinyrna_2023-07-17_20-17-19_summary_stats.csv')
unique(Cel_GA_fc$Classifier)


fc_26g <-Cel_GA_fc[Cel_GA_fc$Classifier %in% c("Genes 26G ", "Pseudogene 26G", "Transposon 26G ", "piRNA 26G", "miRNA 26G"),]
fc_21u<-Cel_GA_fc[Cel_GA_fc$Classifier %in% c("piRNA"),]

rpm <- data.frame(c_26G<- colSums(fc_26g[,c(4:9)]),
                  c_21U<-colSums(fc_21u[,c(4:9)])
                  )
map<-Cel_GA_stat[5,c(2:7)]
rpm_26g<-(rpm$c_26G....colSums.fc_26g...c.4.9.../map)*1000000
rpm_21u<-(rpm$c_21U....colSums.fc_21u...c.4.9.../map)*1000000
rpm2<-melt(rbind(rpm_26g, rpm_21u))
rpm2$type= rep(c("26G", "21U"))
rpm2$stage= rep(c("GA"))
rpm2$cond=rep(c("gtsf-1 (xf251)", "WT"), each= 6)

 a<-rpm2[rpm2$type%in% c("26G"),]
 b<-rpm2[rpm2$type%in% c("21U"),]
 a$cond<-factor(a$cond, levels = c("WT", "gtsf-1 (xf251)"))
 b$cond<-factor(b$cond, levels = c("WT", "gtsf-1 (xf251)"))

 
 #Make RPM plots for Cel_GA
 ggplot(a, aes(x=cond, 
              y=value,
              fill=cond, color = cond)) +
 geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.8)+
#geom_boxplot(width=0.5,outlier.size=0.1) +
  scale_color_manual(values = c("WT" = "steelblue4", 
                                "gtsf-1 (xf251)" =  "darkgreen"
  ))+
  scale_fill_manual(values = c("WT" = "steelblue", 
                               "gtsf-1 (xf251)"=  "forestgreen"
  ))+
  scale_y_continuous(limits = c(0, 5000))+
  ylab(label = '26G RPM') +
  theme_minimal()+
  theme(legend.position = 'top')+
  ggtitle("Levels of 26G RNA reads in C.elegans Gravid Adult samples")
 
 
 
 ##############################
 
 
 
 #For Celegans L4 samples 26G and 21U reads only
 setwd("X:/Shamitha/tinyrna/START_HERE/tinyrna_2023-07-17_21-29-34_Cel_2021_08_L4_26G/tiny-count")
 
 Cel_L4_fc<-read.csv('tinyrna_2023-07-17_21-29-34_feature_counts.csv')
 Cel_L4_stat<-read.csv('tinyrna_2023-07-17_21-29-34_summary_stats.csv')
 unique(Cel_L4_fc$Classifier)
 
 
 fc_26g <-Cel_L4_fc[Cel_L4_fc$Classifier %in% c("Genes 26G ", "Pseudogene 26G", "Transposon 26G ", "piRNA 26G", "miRNA 26G"),]
 fc_21u<-Cel_L4_fc[Cel_L4_fc$Classifier %in% c("piRNA"),]
 
 rpm <- data.frame(c_26G<- colSums(fc_26g[,c(4:9)]),
                   c_21U<-colSums(fc_21u[,c(4:9)])
 )
 map<-Cel_L4_stat[5,c(2:7)]
 rpm_26g<-(rpm$c_26G....colSums.fc_26g...c.4.9.../map)*1000000
 rpm_21u<-(rpm$c_21U....colSums.fc_21u...c.4.9.../map)*1000000
 rpm2<-melt(rbind(rpm_26g, rpm_21u))
 rpm2$type= rep(c("26G", "21U"))
 rpm2$stage= rep(c("L4"))
 rpm2$cond=rep(c("gtsf-1 (xf251)", "WT"), each= 6)
 
 a<-rpm2[rpm2$type%in% c("26G"),]
 b<-rpm2[rpm2$type%in% c("21U"),]
 a$cond<-factor(a$cond, levels = c("WT", "gtsf-1 (xf251)"))
 b$cond<-factor(b$cond, levels = c("WT", "gtsf-1 (xf251)"))
 
 #Make RPM plots for Cel_L4
 ggplot(b, aes(x=cond, 
               y=value,
               fill=cond, color = cond)) +
   geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.8)+
   #geom_boxplot(width=0.5,outlier.size=0.1) +
   scale_color_manual(values = c("WT" = "steelblue4", 
                                 "gtsf-1 (xf251)" =  "darkgreen"
   ))+
   scale_fill_manual(values = c("WT" = "steelblue", 
                                "gtsf-1 (xf251)"=  "forestgreen"
   ))+
   scale_y_continuous(limits = c(95000, 130000))+
   ylab(label = '21U RPM') +
   theme_minimal()+
   theme(legend.position = 'top')+
   ggtitle("Levels of 21U RNA reads in C.elegans L4 samples")
 
 
 ###########
 
 
 
 #briggsae L4 sample 26g 21u
 setwd("X:/Shamitha/tinyrna/START_HERE/tinyrna_2023-07-18_14-02-56_cbr_26G_L4/tiny-count")
 
 Cbr_L4_fc<-read.csv('tinyrna_2023-07-18_14-02-56_feature_counts.csv')
 Cbr_L4_stat<-read.csv('tinyrna_2023-07-18_14-02-56_summary_stats.csv')
 unique(Cbr_L4_fc$Classifier)
 
 
 fc_26g <-Cbr_L4_fc[Cbr_L4_fc$Classifier %in% c("Genes 26G ", "Pseudogene 26G", "Transposon 26G ", "piRNA 26G", "miRNA 26G"),]
 fc_21u<-Cbr_L4_fc[Cbr_L4_fc$Classifier %in% c("piRNA"),]
 
 rpm <- data.frame(c_26G<- colSums(fc_26g[,c(4:9)]),
                   c_21U<-colSums(fc_21u[,c(4:9)])
 )
 map<-Cbr_L4_stat[5,c(2:7)]
 rpm_26g<-(rpm$c_26G....colSums.fc_26g...c.4.9.../map)*1000000
 rpm_21u<-(rpm$c_21U....colSums.fc_21u...c.4.9.../map)*1000000
 rpm2<-melt(rbind(rpm_26g, rpm_21u))
 rpm2$type= rep(c("26G", "21U"))
 rpm2$stage= rep(c("L4"))
 rpm2$cond=rep(c("WT", "gtsf-1 (xf345)"), each= 6)
 
 a<-rpm2[rpm2$type%in% c("26G"),]
 b<-rpm2[rpm2$type%in% c("21U"),]
 a$cond<-factor(a$cond, levels = c("WT", "gtsf-1 (xf345)"))
 b$cond<-factor(b$cond, levels = c("WT", "gtsf-1 (xf345)"))
 
 #Make RPM plots for Cbr_L4
 ggplot(a, aes(x=cond, 
               y=value,
               fill=cond, color = cond)) +
   geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.8)+
   #geom_boxplot(width=0.5,outlier.size=0.1) +
   scale_color_manual(values = c("WT" = "steelblue4", 
                                 "gtsf-1 (xf345)" =  "darkgreen"
   ))+
   scale_fill_manual(values = c("WT" = "steelblue", 
                                "gtsf-1 (xf345)"=  "forestgreen"
   ))+
   scale_y_continuous(limits = c(0, 30000))+
   ylab(label = '26G RPM') +
   theme_minimal()+
   theme(legend.position = 'top')+
   ggtitle("Levels of 26G RNA reads in C.briggsae L4 samples")
 
####################
 
 
 #briggsae Gravid Adult sample 26g 21u
 setwd("X:/Shamitha/tinyrna/START_HERE/tinyrna_2023-07-18_15-19-45_cbr_26G_GA/tiny-count")
 
 Cbr_GA_fc<-read.csv('tinyrna_2023-07-18_15-19-45_feature_counts.csv')
 Cbr_GA_stat<-read.csv('tinyrna_2023-07-18_15-19-45_summary_stats.csv')
 unique(Cbr_GA_fc$Classifier)
 
 
 fc_26g <-Cbr_GA_fc[Cbr_GA_fc$Classifier %in% c("Genes 26G ", "Pseudogene 26G", "Transposon 26G ", "piRNA 26G", "miRNA 26G"),]
 fc_21u<-Cbr_GA_fc[Cbr_GA_fc$Classifier %in% c("piRNA"),]
 
 rpm <- data.frame(c_26G<- colSums(fc_26g[,c(4:9)]),
                   c_21U<-colSums(fc_21u[,c(4:9)])
 )
 map<-Cbr_GA_stat[5,c(2:7)]
 rpm_26g<-(rpm$c_26G....colSums.fc_26g...c.4.9.../map)*1000000
 rpm_21u<-(rpm$c_21U....colSums.fc_21u...c.4.9.../map)*1000000
 rpm2<-melt(rbind(rpm_26g, rpm_21u))
 rpm2$type= rep(c("26G", "21U"))
 rpm2$stage= rep(c("GA"))
 rpm2$cond=rep(c("WT", "gtsf-1 (xf345)"), each= 6)
 
 a<-rpm2[rpm2$type%in% c("26G"),]
 b<-rpm2[rpm2$type%in% c("21U"),]
 a$cond<-factor(a$cond, levels = c("WT", "gtsf-1 (xf345)"))
 b$cond<-factor(b$cond, levels = c("WT", "gtsf-1 (xf345)"))
 
 #Make RPM plots for Cbr_GA
 ggplot(b, aes(x=cond, 
               y=value,
               fill=cond, color = cond)) +
   geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.8)+
   #geom_boxplot(width=0.5,outlier.size=0.1) +
   scale_color_manual(values = c("WT" = "steelblue4", 
                                 "gtsf-1 (xf345)" =  "darkgreen"
   ))+
   scale_fill_manual(values = c("WT" = "steelblue", 
                                "gtsf-1 (xf345)"=  "forestgreen"
   ))+
   scale_y_continuous(limits = c(50000,70000))+
   ylab(label = '21U RPM') +
   theme_minimal()+
   theme(legend.position = 'top')+
   ggtitle("Levels of 21U RNA reads in C.briggsae Gravid Adult samples")
 
 