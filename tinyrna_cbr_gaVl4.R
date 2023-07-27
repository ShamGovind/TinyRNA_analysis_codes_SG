
library(VennDiagram)
setwd("X:/Shamitha/tinyrna/START_HERE/analysis")
targets_l4<-read.csv("26g_targets_cbg_L4.csv")
targets_ga<-read.csv("26g_targets_cbg_GA.csv")

setwd("X:/Shamitha/tinyrna/START_HERE/analysis/briggsae")
gtsf1_targets_l4<-read.csv("gtsf1_targets_cbg_L4.csv")
gtsf1_targets_ga<-read.csv("gtsf1_targets_cbg_GA.csv")



#VennDiagram:
GA_GTSF1 = gtsf1_targets_ga$Feature.ID
GA_26G = targets_ga$Feature.ID

L4_GTSF1 = gtsf1_targets_l4$Feature.ID
L4_26G = targets_l4$Feature.ID

# Create the Venn diagram
venn.diagram(
  x = list(Set1 = GA_GTSF1, Set2 = GA_26G ),
  filename = "venn_cbrigsae_ga_1.png",  # Set filename to NULL for interactive display
  category.names = c("GA_GTSF1", "GA_26G"),
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
  main = "Overlap of 26G RNA targets AND GTSF-1 targets in C. briggsae GA",
  main.cex = 2,
  title = "Sample Title",
  title.col = "black",
  title.cex = 1.8
)
