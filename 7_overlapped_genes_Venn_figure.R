rm(list=ls())
library(RColorBrewer)
library(VennDiagram)
library(scales)
library(gridExtra)  

##################################################BA1########################################################
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/BA1")
enhancer <- read.csv("BA1_annotation.csv")
diff_enhancer <- na.omit(enhancer[[1]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/BA1")
promoter <- read.csv("BA1_annotation.csv")
diff_promoter <- na.omit(promoter[[16]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/RNA/DEG")
RNA <- read.csv("Differential_BA1.csv") 
DEG <- na.omit(RNA[[1]])

setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures/venn")

overlap_1 <- intersect(diff_promoter,diff_enhancer)
overlap_2 <- intersect(diff_promoter, DEG)
overlap_3 <- intersect(diff_enhancer, DEG) 

combined_column <- c(overlap_2, overlap_3)
overlap <- data.frame(Value = combined_column)

write.csv(overlap, file = "BA1_overlap_gene.csv")

figure <- venn.diagram(
  x = list(diff_enhancer, diff_promoter, DEG),
  category.names = c("Enhancer Associated Genes" , "Promoter Associated Genes" , "DEG" ),
  filename = NULL,
  height = 5000 , 
  width = 5000 , 
  resolution = 1000,
  compression = "lzw",
  lwd = 1,
  col= "transparent",
  # col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#8da0cb"), alpha('#66c2a5'), alpha('#fc8d62')),
  cex = 2,
  fontfamily = "sans",
  cat.cex = 1.1,
  # cat.default.pos = "outer",
  # main = c("BA1"), 
  # sub = NULL, main.pos = c(0.5, 1.05), 
  # main.fontface = "bold",
  # main.fontfamily = "sans", main.col = "black",
  # main.cex = 1, main.just = c(0.5, 1), 
  cat.pos = c(-15, 15, 0),
  cat.dist = c(-0.41, -0.41, 0.013),
  cat.fontfamily = "sans",
  cat.col = c("#8da0cb", '#66c2a5', '#fc8d62'),
  rotation.degree = 180)

  pdf(file="venn_BA1.pdf")
  grid.draw(figure)
  dev.off()


##################################################BA2########################################################
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/BA2")
enhancer <- read.csv("BA2_annotation.csv")
diff_enhancer <- na.omit(enhancer[[1]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/BA2")
promoter <- read.csv("BA2_annotation.csv")
diff_promoter <- na.omit(promoter[[16]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/RNA/DEG")
RNA <- read.csv("Differential_BA2.csv") 
DEG <- na.omit(RNA[[1]])

setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures/venn")

overlap_1 <- intersect(diff_promoter,diff_enhancer) 
overlap_2 <- intersect(diff_promoter, DEG) 
overlap_3 <- intersect(diff_enhancer, DEG)  

combined_column <- c(overlap_2, overlap_3)
overlap <- data.frame(Value = combined_column)
write.csv(overlap, file = "BA2_overlap_gene.csv")

figure <- venn.diagram(
  x = list(diff_enhancer, diff_promoter, DEG),
  category.names = c("Enhancer Associated Genes" , "Promoter Associated Genes" , "DEG" ),
  filename = NULL,
  height = 5000 , 
  width = 5000 , 
  resolution = 1000,
  compression = "lzw",
  lwd = 1,
  col= "transparent",
  # col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#8da0cb"), alpha('#66c2a5'), alpha('#fc8d62')),
  cex = 2,
  fontfamily = "sans",
  cat.cex = 1.1,
  # cat.default.pos = "outer",
  # main = c("BA1"), 
  # sub = NULL, main.pos = c(0.5, 1.05), 
  # main.fontface = "bold",
  # main.fontfamily = "sans", main.col = "black",
  # main.cex = 1, main.just = c(0.5, 1), 
  cat.pos = c(-15, 15, 0),
  cat.dist = c(-0.41, -0.41, 0.013),
  cat.fontfamily = "sans",
  cat.col = c("#8da0cb", '#66c2a5', '#fc8d62'),
  rotation.degree = 180)


pdf(file="venn_BA2.pdf")
grid.draw(figure)
dev.off()

##################################################BA3########################################################
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/BA3")
enhancer <- read.csv("BA3_annotation.csv")
diff_enhancer <- na.omit(enhancer[[1]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/BA3")
promoter <- read.csv("BA3_annotation.csv")
diff_promoter <- na.omit(promoter[[16]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/RNA/DEG")
RNA <- read.csv("Differential_BA3.csv") 
DEG <- na.omit(RNA[[1]])

setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures/venn")

overlap_1 <- intersect(diff_promoter,diff_enhancer) 
overlap_2 <- intersect(diff_promoter, DEG) 
overlap_3 <- intersect(diff_enhancer, DEG)  

combined_column <- c(overlap_2, overlap_3)
overlap <- data.frame(Value = combined_column)
write.csv(overlap, file = "BA3_overlap_gene.csv")

figure <- venn.diagram(
  x = list(diff_enhancer, diff_promoter, DEG),
  category.names = c("Enhancer Associated Genes" , "Promoter Associated Genes" , "DEG" ),
  filename = NULL,
  height = 5000 , 
  width = 5000 , 
  resolution = 1000,
  compression = "lzw",
  lwd = 1,
  col= "transparent",
  # col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#8da0cb"), alpha('#66c2a5'), alpha('#fc8d62')),
  cex = 2,
  fontfamily = "sans",
  cat.cex = 1.1,
  # cat.default.pos = "outer",
  # main = c("BA1"), 
  # sub = NULL, main.pos = c(0.5, 1.05), 
  # main.fontface = "bold",
  # main.fontfamily = "sans", main.col = "black",
  # main.cex = 1, main.just = c(0.5, 1), 
  cat.pos = c(-15, 15, 0),
  cat.dist = c(-0.41, -0.41, 0.013),
  cat.fontfamily = "sans",
  cat.col = c("#8da0cb", '#66c2a5', '#fc8d62'),
  rotation.degree = 180)


pdf(file="venn_BA3.pdf")
grid.draw(figure)
dev.off()


##################################################D1C########################################################
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/D1C")
enhancer <- read.csv("D1C_annotation.csv")
diff_enhancer <- na.omit(enhancer[[1]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/D1C")
promoter <- read.csv("D1C_annotation.csv")
diff_promoter <- na.omit(promoter[[16]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/RNA/DEG")
RNA <- read.csv("Differential_D1C.csv") 
DEG <- na.omit(RNA[[1]])

setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures/venn")

overlap_1 <- intersect(diff_promoter,diff_enhancer) 
overlap_2 <- intersect(diff_promoter, DEG) 
overlap_3 <- intersect(diff_enhancer, DEG)  

combined_column <- c(overlap_2, overlap_3)
overlap <- data.frame(Value = combined_column)
write.csv(overlap, file = "D1C_overlap_gene.csv")

figure <- venn.diagram(
  x = list(diff_enhancer, diff_promoter, DEG),
  category.names = c("Enhancer Associated Genes" , "Promoter Associated Genes" , "DEG" ),
  filename = NULL,
  height = 5000 , 
  width = 5000 , 
  resolution = 1000,
  compression = "lzw",
  lwd = 1,
  col= "transparent",
  # col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#8da0cb"), alpha('#66c2a5'), alpha('#fc8d62')),
  cex = 2,
  fontfamily = "sans",
  cat.cex = 1.1,
  # cat.default.pos = "outer",
  # main = c("BA1"), 
  # sub = NULL, main.pos = c(0.5, 1.05), 
  # main.fontface = "bold",
  # main.fontfamily = "sans", main.col = "black",
  # main.cex = 1, main.just = c(0.5, 1), 
  cat.pos = c(-15, 15, 0),
  cat.dist = c(-0.41, -0.41, 0.013),
  cat.fontfamily = "sans",
  cat.col = c("#8da0cb", '#66c2a5', '#fc8d62'),
  rotation.degree = 180)


pdf(file="venn_D1C.pdf")
grid.draw(figure)
dev.off()

##################################################D2C########################################################
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/D2C")
enhancer <- read.csv("D2C_annotation.csv")
diff_enhancer <- na.omit(enhancer[[1]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/D2C")
promoter <- read.csv("D2C_annotation.csv")
diff_promoter <- na.omit(promoter[[16]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/RNA/DEG")
RNA <- read.csv("Differential_D2C.csv") 
DEG <- na.omit(RNA[[1]])

setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures/venn")

overlap_1 <- intersect(diff_promoter,diff_enhancer) 
overlap_2 <- intersect(diff_promoter, DEG) 
overlap_3 <- intersect(diff_enhancer, DEG)  

combined_column <- c(overlap_2, overlap_3)
overlap <- data.frame(Value = combined_column)
write.csv(overlap, file = "D2C_overlap_gene.csv")

figure <- venn.diagram(
  x = list(diff_enhancer, diff_promoter, DEG),
  category.names = c("Enhancer Associated Genes" , "Promoter Associated Genes" , "DEG" ),
  filename = NULL,
  height = 5000 , 
  width = 5000 , 
  resolution = 1000,
  compression = "lzw",
  lwd = 1,
  col= "transparent",
  # col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#8da0cb"), alpha('#66c2a5'), alpha('#fc8d62')),
  cex = 2,
  fontfamily = "sans",
  cat.cex = 1.1,
  # cat.default.pos = "outer",
  # main = c("BA1"), 
  # sub = NULL, main.pos = c(0.5, 1.05), 
  # main.fontface = "bold",
  # main.fontfamily = "sans", main.col = "black",
  # main.cex = 1, main.just = c(0.5, 1), 
  cat.pos = c(-15, 15, 0),
  cat.dist = c(-0.41, -0.41, 0.013),
  cat.fontfamily = "sans",
  cat.col = c("#8da0cb", '#66c2a5', '#fc8d62'),
  rotation.degree = 180)


pdf(file="venn_D2C.pdf")
grid.draw(figure)
dev.off()
##################################################D3C########################################################
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/D3C")
enhancer <- read.csv("D3C_annotation.csv")
diff_enhancer <- na.omit(enhancer[[1]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/D3C")
promoter <- read.csv("D3C_annotation.csv")
diff_promoter <- na.omit(promoter[[16]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/RNA/DEG")
RNA <- read.csv("Differential_D3C.csv") 
DEG <- na.omit(RNA[[1]])

setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures/venn")

overlap_1 <- intersect(diff_promoter,diff_enhancer) 
overlap_2 <- intersect(diff_promoter, DEG) 
overlap_3 <- intersect(diff_enhancer, DEG)  

combined_column <- c(overlap_2, overlap_3)
overlap <- data.frame(Value = combined_column)
write.csv(overlap, file = "D3C_overlap_gene.csv")

figure <- venn.diagram(
  x = list(diff_enhancer, diff_promoter, DEG),
  category.names = c("Enhancer Associated Genes" , "Promoter Associated Genes" , "DEG" ),
  filename = NULL,
  height = 5000 , 
  width = 5000 , 
  resolution = 1000,
  compression = "lzw",
  lwd = 1,
  col= "transparent",
  # col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#8da0cb"), alpha('#66c2a5'), alpha('#fc8d62')),
  cex = 2,
  fontfamily = "sans",
  cat.cex = 1.1,
  # cat.default.pos = "outer",
  # main = c("BA1"), 
  # sub = NULL, main.pos = c(0.5, 1.05), 
  # main.fontface = "bold",
  # main.fontfamily = "sans", main.col = "black",
  # main.cex = 1, main.just = c(0.5, 1), 
  cat.pos = c(-15, 15, 0),
  cat.dist = c(-0.41, -0.41, 0.013),
  cat.fontfamily = "sans",
  cat.col = c("#8da0cb", '#66c2a5', '#fc8d62'),
  rotation.degree = 180)


pdf(file="venn_D3C.pdf")
grid.draw(figure)
dev.off()
##################################################CA1########################################################
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/CA1")
enhancer <- read.csv("CA1_annotation.csv")
diff_enhancer <- na.omit(enhancer[[1]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/CA1")
promoter <- read.csv("CA1_annotation.csv")
diff_promoter <- na.omit(promoter[[16]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/RNA/DEG")
RNA <- read.csv("Differential_CA1.csv") 
DEG <- na.omit(RNA[[1]])

setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures/venn")

overlap_1 <- intersect(diff_promoter,diff_enhancer) 
overlap_2 <- intersect(diff_promoter, DEG) 
overlap_3 <- intersect(diff_enhancer, DEG)  

combined_column <- c(overlap_2, overlap_3)
overlap <- data.frame(Value = combined_column)
write.csv(overlap, file = "CA1_overlap_gene.csv")

figure <- venn.diagram(
  x = list(diff_enhancer, diff_promoter, DEG),
  category.names = c("Enhancer Associated Genes" , "Promoter Associated Genes" , "DEG" ),
  filename = NULL,
  height = 5000 , 
  width = 5000 , 
  resolution = 1000,
  compression = "lzw",
  lwd = 1,
  col= "transparent",
  # col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#8da0cb"), alpha('#66c2a5'), alpha('#fc8d62')),
  cex = 2,
  fontfamily = "sans",
  cat.cex = 1.1,
  # cat.default.pos = "outer",
  # main = c("BA1"), 
  # sub = NULL, main.pos = c(0.5, 1.05), 
  # main.fontface = "bold",
  # main.fontfamily = "sans", main.col = "black",
  # main.cex = 1, main.just = c(0.5, 1), 
  cat.pos = c(-15, 15, 0),
  cat.dist = c(-0.41, -0.41, 0.013),
  cat.fontfamily = "sans",
  cat.col = c("#8da0cb", '#66c2a5', '#fc8d62'),
  rotation.degree = 180)


pdf(file="venn_CA1.pdf")
grid.draw(figure)
dev.off()

##################################################CA2########################################################
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/CA2")
enhancer <- read.csv("CA2_annotation.csv")
diff_enhancer <- na.omit(enhancer[[1]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/CA2")
promoter <- read.csv("CA2_annotation.csv")
diff_promoter <- na.omit(promoter[[16]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/RNA/DEG")
RNA <- read.csv("Differential_CA2.csv") 
DEG <- na.omit(RNA[[1]])

setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures/venn")

overlap_1 <- intersect(diff_promoter,diff_enhancer) 
overlap_2 <- intersect(diff_promoter, DEG) 
overlap_3 <- intersect(diff_enhancer, DEG)  

combined_column <- c(overlap_2, overlap_3)
overlap <- data.frame(Value = combined_column)
write.csv(overlap, file = "CA2_overlap_gene.csv")

figure <- venn.diagram(
  x = list(diff_enhancer, diff_promoter, DEG),
  category.names = c("Enhancer Associated Genes" , "Promoter Associated Genes" , "DEG" ),
  filename = NULL,
  height = 5000 , 
  width = 5000 , 
  resolution = 1000,
  compression = "lzw",
  lwd = 1,
  col= "transparent",
  # col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#8da0cb"), alpha('#66c2a5'), alpha('#fc8d62')),
  cex = 2,
  fontfamily = "sans",
  cat.cex = 1.1,
  # cat.default.pos = "outer",
  # main = c("BA1"), 
  # sub = NULL, main.pos = c(0.5, 1.05), 
  # main.fontface = "bold",
  # main.fontfamily = "sans", main.col = "black",
  # main.cex = 1, main.just = c(0.5, 1), 
  cat.pos = c(-15, 15, 0),
  cat.dist = c(-0.41, -0.41, 0.013),
  cat.fontfamily = "sans",
  cat.col = c("#8da0cb", '#66c2a5', '#fc8d62'),
  rotation.degree = 180)

pdf(file="venn_CA2.pdf")
grid.draw(figure)
dev.off()


##################################################CA3########################################################
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/CA3")
enhancer <- read.csv("CA3_annotation.csv")
diff_enhancer <- na.omit(enhancer[[1]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/CA3")
promoter <- read.csv("CA3_annotation.csv")
diff_promoter <- na.omit(promoter[[16]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/RNA/DEG")
RNA <- read.csv("Differential_CA3.csv") 
DEG <- na.omit(RNA[[1]])

setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures/venn")

overlap_1 <- intersect(diff_promoter,diff_enhancer) 
overlap_2 <- intersect(diff_promoter, DEG) 
overlap_3 <- intersect(diff_enhancer, DEG)  

combined_column <- c(overlap_2, overlap_3)
overlap <- data.frame(Value = combined_column)
write.csv(overlap, file = "CA3_overlap_gene.csv")

figure <- venn.diagram(
  x = list(diff_enhancer, diff_promoter, DEG),
  category.names = c("Enhancer Associated Genes" , "Promoter Associated Genes" , "DEG" ),
  filename = NULL,
  height = 5000 , 
  width = 5000 , 
  resolution = 1000,
  compression = "lzw",
  lwd = 1,
  col= "transparent",
  # col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#8da0cb"), alpha('#66c2a5'), alpha('#fc8d62')),
  cex = 2,
  fontfamily = "sans",
  cat.cex = 1.1,
  # cat.default.pos = "outer",
  # main = c("BA1"), 
  # sub = NULL, main.pos = c(0.5, 1.05), 
  # main.fontface = "bold",
  # main.fontfamily = "sans", main.col = "black",
  # main.cex = 1, main.just = c(0.5, 1), 
  cat.pos = c(-15, 15, 0),
  cat.dist = c(-0.41, -0.41, 0.013),
  cat.fontfamily = "sans",
  cat.col = c("#8da0cb", '#66c2a5', '#fc8d62'),
  rotation.degree = 180)


pdf(file="venn_CA3.pdf")
grid.draw(figure)
dev.off()

##################################################BD1########################################################
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/BD1")
enhancer <- read.csv("BD1_annotation.csv")
diff_enhancer <- na.omit(enhancer[[1]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/BD1")
promoter <- read.csv("BD1_annotation.csv")
diff_promoter <- na.omit(promoter[[16]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/RNA/DEG")
RNA <- read.csv("Differential_BD1.csv") 
DEG <- na.omit(RNA[[1]])

setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures/venn")

overlap_1 <- intersect(diff_promoter,diff_enhancer) 
overlap_2 <- intersect(diff_promoter, DEG) 
overlap_3 <- intersect(diff_enhancer, DEG)  

combined_column <- c(overlap_2, overlap_3)
overlap <- data.frame(Value = combined_column)
write.csv(overlap, file = "BD1_overlap_gene.csv")

figure <- venn.diagram(
  x = list(diff_enhancer, diff_promoter, DEG),
  category.names = c("Enhancer Associated Genes" , "Promoter Associated Genes" , "DEG" ),
  filename = NULL,
  height = 5000 , 
  width = 5000 , 
  resolution = 1000,
  compression = "lzw",
  lwd = 1,
  col= "transparent",
  # col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#8da0cb"), alpha('#66c2a5'), alpha('#fc8d62')),
  cex = 2,
  fontfamily = "sans",
  cat.cex = 1.1,
  # cat.default.pos = "outer",
  # main = c("BA1"), 
  # sub = NULL, main.pos = c(0.5, 1.05), 
  # main.fontface = "bold",
  # main.fontfamily = "sans", main.col = "black",
  # main.cex = 1, main.just = c(0.5, 1), 
  cat.pos = c(-15, 15, 0),
  cat.dist = c(-0.41, -0.41, 0.013),
  cat.fontfamily = "sans",
  cat.col = c("#8da0cb", '#66c2a5', '#fc8d62'),
  rotation.degree = 180)


pdf(file="venn_BD1.pdf")
grid.draw(figure)
dev.off()


##################################################BD2########################################################
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/BD2")
enhancer <- read.csv("BD2_annotation.csv")
diff_enhancer <- na.omit(enhancer[[1]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/BD2")
promoter <- read.csv("BD2_annotation.csv")
diff_promoter <- na.omit(promoter[[16]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/RNA/DEG")
RNA <- read.csv("Differential_BD2.csv") 
DEG <- na.omit(RNA[[1]])

setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures/venn")

overlap_1 <- intersect(diff_promoter,diff_enhancer) 
overlap_2 <- intersect(diff_promoter, DEG) 
overlap_3 <- intersect(diff_enhancer, DEG)  

combined_column <- c(overlap_2, overlap_3)
overlap <- data.frame(Value = combined_column)
write.csv(overlap, file = "BD2_overlap_gene.csv")

figure <- venn.diagram(
  x = list(diff_enhancer, diff_promoter, DEG),
  category.names = c("Enhancer Associated Genes" , "Promoter Associated Genes" , "DEG" ),
  filename = NULL,
  height = 5000 , 
  width = 5000 , 
  resolution = 1000,
  compression = "lzw",
  lwd = 1,
  col= "transparent",
  # col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#8da0cb"), alpha('#66c2a5'), alpha('#fc8d62')),
  cex = 2,
  fontfamily = "sans",
  cat.cex = 1.1,
  # cat.default.pos = "outer",
  # main = c("BA1"), 
  # sub = NULL, main.pos = c(0.5, 1.05), 
  # main.fontface = "bold",
  # main.fontfamily = "sans", main.col = "black",
  # main.cex = 1, main.just = c(0.5, 1), 
  cat.pos = c(-15, 15, 0),
  cat.dist = c(-0.41, -0.41, 0.013),
  cat.fontfamily = "sans",
  cat.col = c("#8da0cb", '#66c2a5', '#fc8d62'),
  rotation.degree = 180)


pdf(file="venn_BD2.pdf")
grid.draw(figure)
dev.off()

##################################################BD3########################################################
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/BD3")
enhancer <- read.csv("BD3_annotation.csv")
diff_enhancer <- na.omit(enhancer[[1]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/BD3")
promoter <- read.csv("BD3_annotation.csv")
diff_promoter <- na.omit(promoter[[16]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/RNA/DEG")
RNA <- read.csv("Differential_BD3.csv") 
DEG <- na.omit(RNA[[1]])

setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures/venn")

overlap_1 <- intersect(diff_promoter,diff_enhancer) 
overlap_2 <- intersect(diff_promoter, DEG) 
overlap_3 <- intersect(diff_enhancer, DEG)  

combined_column <- c(overlap_2, overlap_3)
overlap <- data.frame(Value = combined_column)
write.csv(overlap, file = "BD3_overlap_gene.csv")

figure <- venn.diagram(
  x = list(diff_enhancer, diff_promoter, DEG),
  category.names = c("Enhancer Associated Genes" , "Promoter Associated Genes" , "DEG" ),
  filename = NULL,
  height = 5000 , 
  width = 5000 , 
  resolution = 1000,
  compression = "lzw",
  lwd = 1,
  col= "transparent",
  # col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#8da0cb"), alpha('#66c2a5'), alpha('#fc8d62')),
  cex = 2,
  fontfamily = "sans",
  cat.cex = 1.1,
  # cat.default.pos = "outer",
  # main = c("BA1"), 
  # sub = NULL, main.pos = c(0.5, 1.05), 
  # main.fontface = "bold",
  # main.fontfamily = "sans", main.col = "black",
  # main.cex = 1, main.just = c(0.5, 1), 
  cat.pos = c(-15, 15, 0),
  cat.dist = c(-0.41, -0.41, 0.013),
  cat.fontfamily = "sans",
  cat.col = c("#8da0cb", '#66c2a5', '#fc8d62'),
  rotation.degree = 180)


pdf(file="venn_BD3.pdf")
grid.draw(figure)
dev.off()

##################################################A1E########################################################
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/A1E")
enhancer <- read.csv("A1E_annotation.csv")
diff_enhancer <- na.omit(enhancer[[1]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/A1E")
promoter <- read.csv("A1E_annotation.csv")
diff_promoter <- na.omit(promoter[[16]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/RNA/DEG")
RNA <- read.csv("Differential_A1E.csv") 
DEG <- na.omit(RNA[[1]])

setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures/venn")

overlap_1 <- intersect(diff_promoter,diff_enhancer) 
overlap_2 <- intersect(diff_promoter, DEG) 
overlap_3 <- intersect(diff_enhancer, DEG)  

combined_column <- c(overlap_2, overlap_3)
overlap <- data.frame(Value = combined_column)
write.csv(overlap, file = "A1E_overlap_gene.csv")

figure <- venn.diagram(
  x = list(diff_enhancer, diff_promoter, DEG),
  category.names = c("Enhancer Associated Genes" , "Promoter Associated Genes" , "DEG" ),
  filename = NULL,
  height = 5000 , 
  width = 5000 , 
  resolution = 1000,
  compression = "lzw",
  lwd = 1,
  col= "transparent",
  # col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#8da0cb"), alpha('#66c2a5'), alpha('#fc8d62')),
  cex = 2,
  fontfamily = "sans",
  cat.cex = 1.1,
  # cat.default.pos = "outer",
  # main = c("BA1"), 
  # sub = NULL, main.pos = c(0.5, 1.05), 
  # main.fontface = "bold",
  # main.fontfamily = "sans", main.col = "black",
  # main.cex = 1, main.just = c(0.5, 1), 
  cat.pos = c(-15, 15, 0),
  cat.dist = c(-0.41, -0.41, 0.013),
  cat.fontfamily = "sans",
  cat.col = c("#8da0cb", '#66c2a5', '#fc8d62'),
  rotation.degree = 180)


pdf(file="venn_A1E.pdf")
grid.draw(figure)
dev.off()

##################################################A2E########################################################
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/A2E")
enhancer <- read.csv("A2E_annotation.csv")
diff_enhancer <- na.omit(enhancer[[1]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/A2E")
promoter <- read.csv("A2E_annotation.csv")
diff_promoter <- na.omit(promoter[[16]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/RNA/DEG")
RNA <- read.csv("Differential_A2E.csv") 
DEG <- na.omit(RNA[[1]])

setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures/venn")

overlap_1 <- intersect(diff_promoter,diff_enhancer) 
overlap_2 <- intersect(diff_promoter, DEG) 
overlap_3 <- intersect(diff_enhancer, DEG)  

combined_column <- c(overlap_2, overlap_3)
overlap <- data.frame(Value = combined_column)
write.csv(overlap, file = "A2E_overlap_gene.csv")

figure <- venn.diagram(
  x = list(diff_enhancer, diff_promoter, DEG),
  category.names = c("Enhancer Associated Genes" , "Promoter Associated Genes" , "DEG" ),
  filename = NULL,
  height = 5000 , 
  width = 5000 , 
  resolution = 1000,
  compression = "lzw",
  lwd = 1,
  col= "transparent",
  # col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#8da0cb"), alpha('#66c2a5'), alpha('#fc8d62')),
  cex = 2,
  fontfamily = "sans",
  cat.cex = 1.1,
  # cat.default.pos = "outer",
  # main = c("BA1"), 
  # sub = NULL, main.pos = c(0.5, 1.05), 
  # main.fontface = "bold",
  # main.fontfamily = "sans", main.col = "black",
  # main.cex = 1, main.just = c(0.5, 1), 
  cat.pos = c(-15, 15, 0),
  cat.dist = c(-0.41, -0.41, 0.013),
  cat.fontfamily = "sans",
  cat.col = c("#8da0cb", '#66c2a5', '#fc8d62'),
  rotation.degree = 180)


pdf(file="venn_A2E.pdf")
grid.draw(figure)
dev.off()

##################################################A3E########################################################
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/A3E")
enhancer <- read.csv("A3E_annotation.csv")
diff_enhancer <- na.omit(enhancer[[1]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/A3E")
promoter <- read.csv("A3E_annotation.csv")
diff_promoter <- na.omit(promoter[[16]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/RNA/DEG")
RNA <- read.csv("Differential_A3E.csv") 
DEG <- na.omit(RNA[[1]])

setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures/venn")

overlap_1 <- intersect(diff_promoter,diff_enhancer) 
overlap_2 <- intersect(diff_promoter, DEG) 
overlap_3 <- intersect(diff_enhancer, DEG)  

combined_column <- c(overlap_2, overlap_3)
overlap <- data.frame(Value = combined_column)
write.csv(overlap, file = "A3E_overlap_gene.csv")

figure <- venn.diagram(
  x = list(diff_enhancer, diff_promoter, DEG),
  category.names = c("Enhancer Associated Genes" , "Promoter Associated Genes" , "DEG" ),
  filename = NULL,
  height = 5000 , 
  width = 5000 , 
  resolution = 1000,
  compression = "lzw",
  lwd = 1,
  col= "transparent",
  # col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#8da0cb"), alpha('#66c2a5'), alpha('#fc8d62')),
  cex = 2,
  fontfamily = "sans",
  cat.cex = 1.1,
  # cat.default.pos = "outer",
  # main = c("BA1"), 
  # sub = NULL, main.pos = c(0.5, 1.05), 
  # main.fontface = "bold",
  # main.fontfamily = "sans", main.col = "black",
  # main.cex = 1, main.just = c(0.5, 1), 
  cat.pos = c(-15, 15, 0),
  cat.dist = c(-0.41, -0.41, 0.013),
  cat.fontfamily = "sans",
  cat.col = c("#8da0cb", '#66c2a5', '#fc8d62'),
  rotation.degree = 180)


pdf(file="venn_A3E.pdf")
grid.draw(figure)
dev.off()

##################################################EG########################################################
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/EG")
enhancer <- read.csv("EG_annotation.csv")
diff_enhancer <- na.omit(enhancer[[1]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/EG")
promoter <- read.csv("EG_annotation.csv")
diff_promoter <- na.omit(promoter[[16]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/RNA/DEG")
RNA <- read.csv("Differential_EG.csv") 
DEG <- na.omit(RNA[[1]])

setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures/venn")

overlap_1 <- intersect(diff_promoter,diff_enhancer) 
overlap_2 <- intersect(diff_promoter, DEG) 
overlap_3 <- intersect(diff_enhancer, DEG)  

combined_column <- c(overlap_2, overlap_3)
overlap <- data.frame(Value = combined_column)
write.csv(overlap, file = "EG_overlap_gene.csv")

figure <- venn.diagram(
  x = list(diff_enhancer, diff_promoter, DEG),
  category.names = c("Enhancer Associated Genes" , "Promoter Associated Genes" , "DEG" ),
  filename = NULL,
  height = 5000 , 
  width = 5000 , 
  resolution = 1000,
  compression = "lzw",
  lwd = 1,
  col= "transparent",
  # col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#8da0cb"), alpha('#66c2a5'), alpha('#fc8d62')),
  cex = 2,
  fontfamily = "sans",
  cat.cex = 1.1,
  # cat.default.pos = "outer",
  # main = c("BA1"), 
  # sub = NULL, main.pos = c(0.5, 1.05), 
  # main.fontface = "bold",
  # main.fontfamily = "sans", main.col = "black",
  # main.cex = 1, main.just = c(0.5, 1), 
  cat.pos = c(-15, 15, 0),
  cat.dist = c(-0.41, -0.41, 0.013),
  cat.fontfamily = "sans",
  cat.col = c("#8da0cb", '#66c2a5', '#fc8d62'),
  rotation.degree = 180)


pdf(file="venn_EG.pdf")
grid.draw(figure)
dev.off()

##################################################FE########################################################
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/FE")
enhancer <- read.csv("FE_annotation.csv")
diff_enhancer <- na.omit(enhancer[[1]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/FE")
promoter <- read.csv("FE_annotation.csv")
diff_promoter <- na.omit(promoter[[16]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/RNA/DEG")
RNA <- read.csv("Differential_FE.csv") 
DEG <- na.omit(RNA[[1]])

setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures/venn")

overlap_1 <- intersect(diff_promoter,diff_enhancer) 
overlap_2 <- intersect(diff_promoter, DEG) 
overlap_3 <- intersect(diff_enhancer, DEG)  

combined_column <- c(overlap_2, overlap_3)
overlap <- data.frame(Value = combined_column)
write.csv(overlap, file = "FE_overlap_gene.csv")

figure <- venn.diagram(
  x = list(diff_enhancer, diff_promoter, DEG),
  category.names = c("Enhancer Associated Genes" , "Promoter Associated Genes" , "DEG" ),
  filename = NULL,
  height = 5000 , 
  width = 5000 , 
  resolution = 1000,
  compression = "lzw",
  lwd = 1,
  col= "transparent",
  # col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#8da0cb"), alpha('#66c2a5'), alpha('#fc8d62')),
  cex = 2,
  fontfamily = "sans",
  cat.cex = 1.1,
  # cat.default.pos = "outer",
  # main = c("BA1"), 
  # sub = NULL, main.pos = c(0.5, 1.05), 
  # main.fontface = "bold",
  # main.fontfamily = "sans", main.col = "black",
  # main.cex = 1, main.just = c(0.5, 1), 
  cat.pos = c(-15, 15, 0),
  cat.dist = c(-0.41, -0.41, 0.013),
  cat.fontfamily = "sans",
  cat.col = c("#8da0cb", '#66c2a5', '#fc8d62'),
  rotation.degree = 180)


pdf(file="venn_FE.pdf")
grid.draw(figure)
dev.off()

##################################################FG########################################################
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/FG")
enhancer <- read.csv("FG_annotation.csv")
diff_enhancer <- na.omit(enhancer[[1]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/FG")
promoter <- read.csv("FG_annotation.csv")
diff_promoter <- na.omit(promoter[[16]])
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/RNA/DEG")
RNA <- read.csv("Differential_FG.csv") 
DEG <- na.omit(RNA[[1]])

setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures/venn")

overlap_1 <- intersect(diff_promoter,diff_enhancer) 
overlap_2 <- intersect(diff_promoter, DEG) 
overlap_3 <- intersect(diff_enhancer, DEG)  

combined_column <- c(overlap_2, overlap_3)
overlap <- data.frame(Value = combined_column)

write.csv(overlap, file = "FG_overlap_gene.csv")

figure <- venn.diagram(
  x = list(diff_enhancer, diff_promoter, DEG),
  category.names = c("Enhancer Associated Genes" , "Promoter Associated Genes" , "DEG" ),
  filename = NULL,
  height = 5000 , 
  width = 5000 , 
  resolution = 1000,
  compression = "lzw",
  lwd = 1,
  col= "transparent",
  # col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#8da0cb"), alpha('#66c2a5'), alpha('#fc8d62')),
  cex = 2,
  fontfamily = "sans",
  cat.cex = 1.1,
  # cat.default.pos = "outer",
  # main = c("BA1"), 
  # sub = NULL, main.pos = c(0.5, 1.05), 
  # main.fontface = "bold",
  # main.fontfamily = "sans", main.col = "black",
  # main.cex = 1, main.just = c(0.5, 1), 
  cat.pos = c(-15, 15, 0),
  cat.dist = c(-0.41, -0.41, 0.013),
  cat.fontfamily = "sans",
  cat.col = c("#8da0cb", '#66c2a5', '#fc8d62'),
  rotation.degree = 180)


pdf(file="venn_FG.pdf")
grid.draw(figure)
dev.off()
