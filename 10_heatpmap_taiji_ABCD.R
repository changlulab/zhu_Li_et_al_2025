###########prep#############
rm(list = ls())
options(stringsAsFactors = F)

# Load required libraries
library(dplyr)
library(MASS)
library(fitdistrplus)
library(reshape2)
library(ggplot2)
library(DESeq2)
library(matrixStats)
library(ggpubr)
library(ggh4x)
library(viridis)

# Set working directory
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/taiji/output_from_taiji/")

# Read in data
taiji_res = read.csv('GeneRanks_ABCD.tsv', sep = '\t', header = T)
taiji_res_E = read.csv('GeneRanks_EFG.tsv', sep = '\t', header = T)

###################################BA1######
# Select relevant columns and rename them
new_BA1 = bind_cols(
  taiji_res[,30], taiji_res[,32], taiji_res[,34], taiji_res[,36], taiji_res[,38], taiji_res[,40],
  taiji_res[,2], taiji_res[,4], taiji_res[,6], taiji_res[,7], taiji_res[,9], taiji_res[,11]
)
colnames(new_BA1) = c(
  'GroupB_1', 'GroupB_2', 'GroupB_3', 'GroupB_4', 'GroupB_5', 'GroupB_6',
  'GroupA_1', 'GroupA_2', 'GroupA_3', 'GroupA_4', 'GroupA_5', 'GroupA_6'
)
new_BA1$TFs = taiji_res$X

# Z-score normalization
zscore_BA1 = new_BA1
for (i in 1:(ncol(new_BA1) - 1)) {
  zscore_BA1[, i] = scale(new_BA1[, i])
}
zscore_BA1$mean = rowMeans(zscore_BA1[, 1:12])
zscore_BA1$sd = rowSds(as.matrix(zscore_BA1[, 1:12]))

# Standardize again by row mean and sd
tmp = data.frame((as.matrix(zscore_BA1[,1:12]) - zscore_BA1$mean) / zscore_BA1$sd)
tmp$TFs = taiji_res$X

# Read differential TFs
diff_TF = read.csv('Diff_TFs_BA1.csv')
increase_TFs = diff_TF[which(diff_TF$diff > 0), ]
decreade_TFs = diff_TF[which(diff_TF$diff < 0), ]

# Filter matrix
TFs_matrix_up = tmp[which(tmp$TFs %in% increase_TFs$TFs), ]
TFs_matrix_down = tmp[which(tmp$TFs %in% decreade_TFs$TFs), ]
TFs_matrix = bind_rows(TFs_matrix_up, TFs_matrix_down)

# Prepare data for heatmap
new = as.matrix(TFs_matrix[, 1:12])
rownames(new) = TFs_matrix$TFs
melted_BA1 = melt(new)

# Annotate groups
melted_BA1$Group <- NA
melted_BA1$Group[1:54] <- "GroupB"
melted_BA1$Group[55:108] <- "GroupA1"

# Plot heatmap
pdf(file = "heatmap_BA1_key.pdf", width = 5, height = 3)
ggplot(data = melted_BA1, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() + 
  scale_fill_viridis(name = "PRP z-score") +
  facet_grid2(. ~ Group, scales = "free_x", space = "free_x", switch = "x") +
  guides(fill = guide_colourbar(title = "PRP z-score"), color = "black") +
  theme(
    axis.title.x     = element_blank(),
    axis.title.y     = element_blank(),
    axis.text.x      = element_blank(),
    panel.spacing.x  = unit(0, "mm"),
    panel.spacing.y  = unit(0, "mm"),
    panel.border     = element_rect(fill = NA, color = "black", size = 1),
    panel.background = element_blank(),
    panel.grid       = element_blank(),
    strip.background = element_rect(fill = NA),
    strip.placement  = "outside",
    legend.background = element_rect(fill = NA, size = 0.5),
    legend.key       = element_rect(colour = "black", size = 0.5),
    plot.margin      = unit(c(0.5, 1, 0.5, 1), "cm"),
    legend.position  = "bottom"
  )
dev.off()
###################################BA2######
new_BA2 = bind_cols(taiji_res[,30],taiji_res[,32],taiji_res[,34],taiji_res[,36],taiji_res[,38],taiji_res[,40],
                    taiji_res[,13],taiji_res[,15],taiji_res[,17])
colnames(new_BA2) = c('GroupB_1','GroupB_2','GroupB_3','GroupB_4',
                      'GroupB_5','GroupB_6',
                      'GroupA_7','GroupA_8','GroupA_9')
new_BA2$TFs =taiji_res$X
zscore_BA2 = new_BA2
for (i in 1:(ncol(new_BA2)-1)) {
  zscore_BA2[,i] = scale(new_BA2[,i])
}
zscore_BA2$mean = rowMeans(zscore_BA2[,1:9])
zscore_BA2$sd = rowSds(as.matrix(zscore_BA2[,1:9]))
tmp = data.frame((as.matrix(zscore_BA2[,1:9])-zscore_BA2$mean)/zscore_BA2$sd)
tmp$TFs = taiji_res$X
# setwd("~/Documents/Virginia Tech/Project_3_MIA mouse/taiji/ABCD/BA/")
# TFs_B = read.csv('potential_key_TFs_BA2_GroupB2.csv') 
# TFs_A = read.csv('potential_key_TFs_BA2_GroupA2.csv')
# TFs = union(TFs_B$taiji_res.X,TFs_A$taiji_res.X)
# setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/taiji/output_from_taiji/")
diff_TF = read.csv('Diff_TFs_BA2.csv')
# diff_TF = diff_TF[which(diff_TF$TFs %in% TFs),]
increase_TFs = diff_TF[which(diff_TF$diff > 0),]
decreade_TFs =diff_TF[which(diff_TF$diff < 0),]
# up_TF = TFs[which(TFs$V1 %in% increase_TFs$TFs),]
# down_TF = TFs[which(TFs$V1 %in% decreade_TFs$TFs),]
# write.table(up_TF,'up_TFs.txt',row.names = FALSE,quote = FALSE)
# write.table(down_TF,'down_TF.txt',row.names = FALSE,quote = FALSE)
TFs_matrix_up = tmp[which(tmp$TFs %in% increase_TFs$TFs),]
TFs_matrix_down = tmp[which(tmp$TFs %in% decreade_TFs$TFs),]
TFs_matrix = bind_rows(TFs_matrix_up,TFs_matrix_down)
new = as.matrix(TFs_matrix[,1:9])
rownames(new)= TFs_matrix$TFs
melted_BA2 = melt(new)

melted_BA2$Group <- NA
melted_BA2$Group[1:36] <- "GroupB"
melted_BA2$Group[37:54] <- "GroupA2"
pdf(file = "heatmap_BA2_key.pdf",width = 5,height = 4)
ggplot(data = melted_BA2, aes(x=Var2, y=Var1, fill= value)) +
  geom_tile() + 
  scale_fill_viridis(name = "PRP z-score") +
  facet_grid2(. ~ Group, scales = "free_x", space = "free_x", switch = "x") +
  guides(fill = guide_colourbar(title = "PRP z-score"), color = "black") +
  theme(
    axis.title.x     = element_blank(),
    axis.title.y     = element_blank(),
    axis.text.x      = element_blank(),
    panel.spacing.x  = unit(0, "mm"),
    panel.spacing.y  = unit(0, "mm"),
    panel.border     = element_rect(fill = NA, color = "black", size = 1),
    panel.background = element_blank(),
    panel.grid       = element_blank(),
    strip.background = element_rect(fill = NA),
    strip.placement  = "outside",
    legend.background = element_rect(fill = NA, size = 0.5),
    legend.key       = element_rect(colour = "black", size = 0.5),
    plot.margin      = unit(c(0.5, 1, 0.5, 1), "cm"),
    legend.position  = "bottom"
  )
dev.off()
###################################BA3######
new_BA3 = bind_cols(taiji_res[,30],taiji_res[,32],taiji_res[,34],taiji_res[,36],taiji_res[,38],taiji_res[,40],
                    taiji_res[,19],taiji_res[,20],taiji_res[,22],taiji_res[,24],
                    taiji_res[,26],taiji_res[,28])
colnames(new_BA3) = c('GroupB_1','GroupB_2','GroupB_3','GroupB_4',
                      'GroupB_5','GroupB_6',
                      'GroupA_10','GroupA_11','GroupA_12','GroupA_13',
                      'GroupA_14','GroupA_15')
new_BA3$TFs =taiji_res$X
zscore_BA3 = new_BA3
for (i in 1:(ncol(new_BA3)-1)) {
  zscore_BA3[,i] = scale(new_BA3[,i])
}
zscore_BA3$mean = rowMeans(zscore_BA3[,1:12])
zscore_BA3$sd = rowSds(as.matrix(zscore_BA3[,1:12]))
tmp = data.frame((as.matrix(zscore_BA3[,1:12])-zscore_BA3$mean)/zscore_BA3$sd)
tmp$TFs = taiji_res$X
# setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/taiji/output_from_taiji/")
# TFs_B = read.csv('potential_key_TFs_BA3_GroupB3.csv') 
# TFs_A = read.csv('potential_key_TFs_BA3_GroupA3.csv')
# TFs = union(TFs_B$taiji_res.X,TFs_A$taiji_res.X)
# setwd("~/Documents/Virginia Tech/Project_3_MIA mouse/taiji/ABCD/BA/")
diff_TF = read.csv('Diff_TFs_BA3.csv')
# diff_TF = diff_TF[which(diff_TF$TFs %in% TFs),]
increase_TFs = diff_TF[which(diff_TF$diff > 0),]
decreade_TFs =diff_TF[which(diff_TF$diff < 0),]
# up_TF = TFs[which(TFs$V1 %in% increase_TFs$TFs),]
# down_TF = TFs[which(TFs$V1 %in% decreade_TFs$TFs),]
# write.table(up_TF,'up_TFs.txt',row.names = FALSE,quote = FALSE)
# write.table(down_TF,'down_TF.txt',row.names = FALSE,quote = FALSE)
TFs_matrix_up = tmp[which(tmp$TFs %in% increase_TFs$TFs),]
TFs_matrix_down = tmp[which(tmp$TFs %in% decreade_TFs$TFs),]
TFs_matrix = bind_rows(TFs_matrix_up,TFs_matrix_down)
new = as.matrix(TFs_matrix[,1:12])
rownames(new)= TFs_matrix$TFs
melted_BA3 = melt(new)

melted_BA3$Group <- NA
melted_BA3$Group[1:42] <- "GroupB"
melted_BA3$Group[43:84] <- "GroupA3"
pdf(file = "heatmap_BA3_key.pdf",width = 5,height = 4)
ggplot(data = melted_BA3, aes(x=Var2, y=Var1, fill= value)) +
  geom_tile() + 
  scale_fill_viridis(name = "PRP z-score") +
  facet_grid2(. ~ Group, scales = "free_x", space = "free_x", switch = "x") +
  guides(fill = guide_colourbar(title = "PRP z-score"), color = "black") +
  theme(
    axis.title.x     = element_blank(),
    axis.title.y     = element_blank(),
    axis.text.x      = element_blank(),
    panel.spacing.x  = unit(0, "mm"),
    panel.spacing.y  = unit(0, "mm"),
    panel.border     = element_rect(fill = NA, color = "black", size = 1),
    panel.background = element_blank(),
    panel.grid       = element_blank(),
    strip.background = element_rect(fill = NA),
    strip.placement  = "outside",
    legend.background = element_rect(fill = NA, size = 0.5),
    legend.key       = element_rect(colour = "black", size = 0.5),
    plot.margin      = unit(c(0.5, 1, 0.5, 1), "cm"),
    legend.position  = "bottom"
  )
dev.off()
###################################D1C######
new_D1C = bind_cols(taiji_res[,54],taiji_res[,56],taiji_res[,58],taiji_res[,60],taiji_res[,62],taiji_res[,64],
                    taiji_res[,42],taiji_res[,44],taiji_res[,46],taiji_res[,48],taiji_res[,50],taiji_res[,52])
colnames(new_D1C) = c('GroupD_1','GroupD_2','GroupD_3','GroupD_4','GroupD_5','GroupD_6',
                      'GroupC_1','GroupC_2','GroupC_3','GroupC_4','GroupC_5','GroupC_6')
new_D1C$TFs =taiji_res$X
zscore_D1C = new_D1C
for (i in 1:(ncol(new_D1C)-1)) {
  zscore_D1C[,i] = scale(new_D1C[,i])
}
zscore_D1C$mean = rowMeans(zscore_D1C[,1:12])
zscore_D1C$sd = rowSds(as.matrix(zscore_D1C[,1:12]))
tmp = data.frame((as.matrix(zscore_D1C[,1:12])-zscore_D1C$mean)/zscore_D1C$sd)
tmp$TFs = taiji_res$X
diff_TF = read.csv('Diff_TFs_D1C.csv')

increase_TFs = diff_TF[which(diff_TF$diff > 0),]
decreade_TFs =diff_TF[which(diff_TF$diff < 0),]

TFs_matrix_up = tmp[which(tmp$TFs %in% increase_TFs$TFs),]
TFs_matrix_down = tmp[which(tmp$TFs %in% decreade_TFs$TFs),]
TFs_matrix = bind_rows(TFs_matrix_up,TFs_matrix_down)
new = as.matrix(TFs_matrix[,1:12])
rownames(new)= TFs_matrix$TFs
melted_D1C = melt(new)

melted_D1C$Group <- NA
melted_D1C$Group[1:78] <- "GroupD1"
melted_D1C$Group[79:156] <- "GroupC"
pdf(file = "heatmap_D1C_key.pdf",width = 5,height = 3)
ggplot(data = melted_D1C, aes(x=Var2, y=Var1, fill= value)) +
  geom_tile() + scale_fill_viridis(name = "PRP z-score") +
  facet_grid2(.~Group, scales = "free_x", space = "free_x", switch = "x")+
  guides(fill = guide_colourbar(title = "PRP z-score"), color= "black")+
  theme(axis.title.x     = element_blank(),
        axis.title.y     = element_blank(),
        axis.text.x      = element_blank(),
        panel.spacing.x  = unit(0, "mm"),
        panel.spacing.y  = unit(0, "mm"),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background = element_blank(),
        panel.grid       = element_blank(),
        strip.background = element_rect(fill= NA),
        strip.placement = "outside",
        legend.background = element_rect(fill=NA,
                                         size=0.5),
        legend.key = element_rect(colour = "black", size=0.5),
        plot.margin = unit(c(0.5,1,0.5,1), "cm"),
        legend.position = "bottom")

dev.off()
###################################D2C######
new_D2C = bind_cols(taiji_res[,66],taiji_res[,68],taiji_res[,70],
                    taiji_res[,42],taiji_res[,44],taiji_res[,46],taiji_res[,48],taiji_res[,50],taiji_res[,52])
colnames(new_D2C) = c('GroupD_7','GroupD_8','GroupD_9',
                      'GroupC_1','GroupC_2','GroupC_3','GroupC_4','GroupC_5','GroupC_6')
new_D2C$TFs =taiji_res$X
zscore_D2C = new_D2C
for (i in 1:(ncol(new_D2C)-1)) {
  zscore_D2C[,i] = scale(new_D2C[,i])
}
zscore_D2C$mean = rowMeans(zscore_D2C[,1:9])
zscore_D2C$sd = rowSds(as.matrix(zscore_D2C[,1:9]))
tmp = data.frame((as.matrix(zscore_D2C[,1:9])-zscore_D2C$mean)/zscore_D2C$sd)
tmp$TFs = taiji_res$X
diff_TF = read.csv('Diff_TFs_D2C.csv')

increase_TFs = diff_TF[which(diff_TF$diff > 0),]
decreade_TFs =diff_TF[which(diff_TF$diff < 0),]

TFs_matrix_up = tmp[which(tmp$TFs %in% increase_TFs$TFs),]
TFs_matrix_down = tmp[which(tmp$TFs %in% decreade_TFs$TFs),]
TFs_matrix = bind_rows(TFs_matrix_up,TFs_matrix_down)
new = as.matrix(TFs_matrix[,1:9])
rownames(new)= TFs_matrix$TFs
melted_D2C = melt(new)

melted_D2C$Group <- NA
melted_D2C$Group[1:21] <- "GroupD2"
melted_D2C$Group[22:63] <- "GroupC"
pdf(file = "heatmap_D2C_key.pdf",width = 5,height = 3)
ggplot(data = melted_D2C, aes(x=Var2, y=Var1, fill= value)) +
  geom_tile() + scale_fill_viridis(name = "PRP z-score") +
  facet_grid2(.~Group, scales = "free_x", space = "free_x", switch = "x")+
  guides(fill = guide_colourbar(title = "PRP z-score"), color= "black")+
  theme(axis.title.x     = element_blank(),
        axis.title.y     = element_blank(),
        axis.text.x      = element_blank(),
        panel.spacing.x  = unit(0, "mm"),
        panel.spacing.y  = unit(0, "mm"),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background = element_blank(),
        panel.grid       = element_blank(),
        strip.background = element_rect(fill= NA),
        strip.placement = "outside",
        legend.background = element_rect(fill=NA,
                                         size=0.5),
        legend.key = element_rect(colour = "black", size=0.5),
        plot.margin = unit(c(0.5,1,0.5,1), "cm"),
        legend.position = "bottom")

dev.off()
###################################D3C######
new_D3C = bind_cols(taiji_res[,72],taiji_res[,74],taiji_res[,76],taiji_res[,77],taiji_res[,79],
                    taiji_res[,42],taiji_res[,44],taiji_res[,46],taiji_res[,48],taiji_res[,50],taiji_res[,52])
colnames(new_D3C) = c('GroupD_11','GroupD_12','GroupD_13','GroupD_14','GroupD_15',
                      'GroupC_1','GroupC_2','GroupC_3','GroupC_4','GroupC_5','GroupC_6')
new_D3C$TFs =taiji_res$X
zscore_D3C = new_D3C
for (i in 1:(ncol(new_D3C)-1)) {
  zscore_D3C[,i] = scale(new_D3C[,i])
}
zscore_D3C$mean = rowMeans(zscore_D3C[,1:11])
zscore_D3C$sd = rowSds(as.matrix(zscore_D3C[,1:11]))
tmp = data.frame((as.matrix(zscore_D3C[,1:11])-zscore_D3C$mean)/zscore_D3C$sd)
tmp$TFs = taiji_res$X
diff_TF = read.csv('Diff_TFs_D3C.csv')

increase_TFs = diff_TF[which(diff_TF$diff > 0),]
decreade_TFs =diff_TF[which(diff_TF$diff < 0),]

TFs_matrix_up = tmp[which(tmp$TFs %in% increase_TFs$TFs),]
TFs_matrix_down = tmp[which(tmp$TFs %in% decreade_TFs$TFs),]
TFs_matrix = bind_rows(TFs_matrix_up,TFs_matrix_down)
new = as.matrix(TFs_matrix[,1:11])
rownames(new)= TFs_matrix$TFs
melted_D3C = melt(new)

melted_D3C$Group <- NA
melted_D3C$Group[1:60] <- "GroupD3"
melted_D3C$Group[61:132] <- "GroupC"
pdf(file = "heatmap_D3C_key.pdf",width = 5,height = 3)
ggplot(data = melted_D3C, aes(x=Var2, y=Var1, fill= value)) +
  geom_tile() + scale_fill_viridis(name = "PRP z-score") +
  facet_grid2(.~Group, scales = "free_x", space = "free_x", switch = "x")+
  guides(fill = guide_colourbar(title = "PRP z-score"), color= "black")+
  theme(axis.title.x     = element_blank(),
        axis.title.y     = element_blank(),
        axis.text.x      = element_blank(),
        panel.spacing.x  = unit(0, "mm"),
        panel.spacing.y  = unit(0, "mm"),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background = element_blank(),
        panel.grid       = element_blank(),
        strip.background = element_rect(fill= NA),
        strip.placement = "outside",
        legend.background = element_rect(fill=NA,
                                         size=0.5),
        legend.key = element_rect(colour = "black", size=0.5),
        plot.margin = unit(c(0.5,1,0.5,1), "cm"),
        legend.position = "bottom")

dev.off()
###################################CA1######
new_CA1 = bind_cols(taiji_res[,42],taiji_res[,44],taiji_res[,46],taiji_res[,48],taiji_res[,50],taiji_res[,52],
                    taiji_res[,2],taiji_res[,4],taiji_res[,6],taiji_res[,7],taiji_res[,9],taiji_res[,11])
colnames(new_CA1) = c('GroupC_1','GroupC_2','GroupC_3','GroupC_4','GroupC_5','GroupC_6',
                      'GroupA_1','GroupA_2','GroupA_3','GroupA_4','GroupA_5','GroupA_6')
new_CA1$TFs =taiji_res$X
zscore_CA1 = new_CA1
for (i in 1:(ncol(new_CA1)-1)) {
  zscore_CA1[,i] = scale(new_CA1[,i])
}
zscore_CA1$mean = rowMeans(zscore_CA1[,1:12])
zscore_CA1$sd = rowSds(as.matrix(zscore_CA1[,1:12]))
tmp = data.frame((as.matrix(zscore_CA1[,1:12])-zscore_CA1$mean)/zscore_CA1$sd)
tmp$TFs = taiji_res$X
# setwd("~/Documents/Virginia Tech/Project_3_MIA mouse/taiji/ABCD/CA/")
# TFs_C = read.csv('potential_key_TFs_CA1_GroupC.csv') 
# TFs_A = read.csv('potential_key_TFs_CA1_GroupA1.csv')
# TFs = union(TFs_C$taiji_res.X,TFs_A$taiji_res.X)
# setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/taiji/output_from_taiji/")
diff_TF = read.csv('Diff_TFs_CA1.csv')
# diff_TF = diff_TF[which(diff_TF$TFs %in% TFs),]
increase_TFs = diff_TF[which(diff_TF$diff > 0),]
decreade_TFs =diff_TF[which(diff_TF$diff < 0),]
# up_TF = TFs[which(TFs$V1 %in% increase_TFs$TFs),]
# down_TF = TFs[which(TFs$V1 %in% decreade_TFs$TFs),]
# write.table(up_TF,'up_TFs.txt',row.names = FALSE,quote = FALSE)
# write.table(down_TF,'down_TF.txt',row.names = FALSE,quote = FALSE)
TFs_matrix_up = tmp[which(tmp$TFs %in% increase_TFs$TFs),]
TFs_matrix_down = tmp[which(tmp$TFs %in% decreade_TFs$TFs),]
TFs_matrix = bind_rows(TFs_matrix_up,TFs_matrix_down)
new = as.matrix(TFs_matrix[,1:12])
rownames(new)= TFs_matrix$TFs
melted_CA1 = melt(new)
melted_CA1$Group <- NA
melted_CA1$Group[1:48] <- "GroupC"
melted_CA1$Group[49:96] <- "GroupA1"
pdf(file = "heatmap_CA1_key.pdf",width = 5,height = 3)
ggplot(data = melted_CA1, aes(x=Var2, y=Var1, fill= value)) + 
  geom_tile() + scale_fill_viridis(name = "PRP z-score") +
  facet_grid2(.~Group, scales = "free_x", space = "free_x", switch = "x")+
  guides(fill = guide_colourbar(title = "PRP z-score"), color= "black")+
  theme(axis.title.x     = element_blank(),
        axis.title.y     = element_blank(),
        axis.text.x      = element_blank(),
        panel.spacing.x  = unit(0, "mm"),
        panel.spacing.y  = unit(0, "mm"),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background = element_blank(),
        panel.grid       = element_blank(),
        strip.background = element_rect(fill= NA),
        strip.placement = "outside",
        legend.background = element_rect(fill=NA,
                                         size=0.5),
        plot.margin = unit(c(0.5,1,0.5,1), "cm"),
        legend.position = "bottom")

        

dev.off()


###################################CA2######
new_CA2 = bind_cols(taiji_res[,42],taiji_res[,44],taiji_res[,46],taiji_res[,48],taiji_res[,50],taiji_res[,52],
                    taiji_res[,13],taiji_res[,15],taiji_res[,17])
colnames(new_CA2) = c('GroupC_1','GroupC_2','GroupC_3','GroupC_4','GroupC_5','GroupC_6',
                      'GroupA_7','GroupA_8','GroupA_9')
new_CA2$TFs =taiji_res$X
zscore_CA2 = new_CA2
for (i in 1:(ncol(new_CA2)-1)) {
  zscore_CA2[,i] = scale(new_CA2[,i])
}
zscore_CA2$mean = rowMeans(zscore_CA2[,1:9])
zscore_CA2$sd = rowSds(as.matrix(zscore_CA2[,1:9]))
tmp = data.frame((as.matrix(zscore_CA2[,1:9])-zscore_CA2$mean)/zscore_CA2$sd)
tmp$TFs = taiji_res$X
# setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/taiji/output_from_taiji/")
# TFs_C = read.csv('potential_key_TFs_CA2_GroupC.csv') 
# TFs_A = read.csv('potential_key_TFs_CA2_GroupA2.csv')
# TFs = union(TFs_C$taiji_res.X,TFs_A$taiji_res.X)
# setwd("~/Documents/Virginia Tech/Project_3_MIA mouse/taiji/ABCD/CA/")
diff_TF = read.csv('Diff_TFs_CA2.csv')
# diff_TF = diff_TF[which(diff_TF$TFs %in% TFs),]
increase_TFs = diff_TF[which(diff_TF$diff > 0),]
decreade_TFs =diff_TF[which(diff_TF$diff < 0),]
# up_TF = TFs[which(TFs$V1 %in% increase_TFs$TFs),]
# down_TF = TFs[which(TFs$V1 %in% decreade_TFs$TFs),]
# write.table(up_TF,'up_TFs.txt',row.names = FALSE,quote = FALSE)
# write.table(down_TF,'down_TF.txt',row.names = FALSE,quote = FALSE)
TFs_matrix_up = tmp[which(tmp$TFs %in% increase_TFs$TFs),]
TFs_matrix_down = tmp[which(tmp$TFs %in% decreade_TFs$TFs),]
TFs_matrix = bind_rows(TFs_matrix_up,TFs_matrix_down)
new = as.matrix(TFs_matrix[,1:9])
rownames(new)= TFs_matrix$TFs
melted_CA2 = melt(new)

melted_CA2$Group <- NA
melted_CA2$Group[1:42] <- "GroupC"
melted_CA2$Group[43:63] <- "GroupA2"
pdf(file = "heatmap_CA2_key.pdf",width = 5,height = 3)
ggplot(data = melted_CA2, aes(x=Var2, y=Var1, fill= value)) + 
  geom_tile() + scale_fill_viridis(name = "PRP z-score") +
  facet_grid2(.~Group, scales = "free_x", space = "free_x", switch = "x")+
  guides(fill = guide_colourbar(title = "PRP z-score"), color= "black")+
  theme(axis.title.x     = element_blank(),
        axis.title.y     = element_blank(),
        axis.text.x      = element_blank(),
        panel.spacing.x  = unit(0, "mm"),
        panel.spacing.y  = unit(0, "mm"),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background = element_blank(),
        panel.grid       = element_blank(),
        strip.background = element_rect(fill= NA),
        strip.placement = "outside",
        legend.background = element_rect(fill=NA,
                                         size=0.5),
        plot.margin = unit(c(0.5,1,0.5,1), "cm"),
        legend.position = "bottom")



dev.off()

###################################CA3######
new_CA3 = bind_cols(taiji_res[,42],taiji_res[,44],taiji_res[,46],taiji_res[,48],taiji_res[,50],taiji_res[,52],
                    taiji_res[,19],taiji_res[,20],taiji_res[,22],taiji_res[,24],taiji_res[,26],taiji_res[,28])
colnames(new_CA3) = c('GroupC_1','GroupC_2','GroupC_3','GroupC_4','GroupC_5','GroupC_6',
                      'GroupA_10','GroupA_11','GroupA_12','GroupA_13','GroupA_14','GroupA_15')
new_CA3$TFs =taiji_res$X
zscore_CA3 = new_CA3
for (i in 1:(ncol(new_CA3)-1)) {
  zscore_CA3[,i] = scale(new_CA3[,i])
}
zscore_CA3$mean = rowMeans(zscore_CA3[,1:12])
zscore_CA3$sd = rowSds(as.matrix(zscore_CA3[,1:12]))
tmp = data.frame((as.matrix(zscore_CA3[,1:12])-zscore_CA3$mean)/zscore_CA3$sd)
tmp$TFs = taiji_res$X
# setwd("~/Documents/Virginia Tech/Project_3_MIA mouse/taiji_new/ABCD/CA3")
# TFs_C = read.csv('potential_key_TFs_CA3_GroupC.csv') 
# TFs_A = read.csv('potential_key_TFs_CA3_GroupA3.csv')
# TFs = union(TFs_C$taiji_res.X,TFs_A$taiji_res.X)
# setwd("~/Documents/Virginia Tech/Project_3_MIA mouse/taiji/ABCD/CA/")
diff_TF = read.csv('Diff_TFs_CA3.csv')
# diff_TF = diff_TF[which(diff_TF$TFs %in% TFs),]
increase_TFs = diff_TF[which(diff_TF$diff > 0),]
decreade_TFs =diff_TF[which(diff_TF$diff < 0),]
# up_TF = TFs[which(TFs$V1 %in% increase_TFs$TFs),]
# down_TF = TFs[which(TFs$V1 %in% decreade_TFs$TFs),]
# write.table(up_TF,'up_TFs.txt',row.names = FALSE,quote = FALSE)
# write.table(down_TF,'down_TF.txt',row.names = FALSE,quote = FALSE)
TFs_matrix_up = tmp[which(tmp$TFs %in% increase_TFs$TFs),]
TFs_matrix_down = tmp[which(tmp$TFs %in% decreade_TFs$TFs),]
TFs_matrix = bind_rows(TFs_matrix_up,TFs_matrix_down)
new = as.matrix(TFs_matrix[,1:12])
rownames(new)= TFs_matrix$TFs
melted_CA3 = melt(new)

melted_CA3$Group <- NA
melted_CA3$Group[1:54] <- "GroupC"
melted_CA3$Group[55:108] <- "GroupA3"
pdf(file = "heatmap_CA3_key.pdf",width = 5,height = 3)
ggplot(data = melted_CA3, aes(x=Var2, y=Var1, fill= value)) + 
  geom_tile() + scale_fill_viridis(name = "PRP z-score") +
  facet_grid2(.~Group, scales = "free_x", space = "free_x", switch = "x")+
  guides(fill = guide_colourbar(title = "PRP z-score"), color= "black")+
  theme(axis.title.x     = element_blank(),
        axis.title.y     = element_blank(),
        axis.text.x      = element_blank(),
        panel.spacing.x  = unit(0, "mm"),
        panel.spacing.y  = unit(0, "mm"),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background = element_blank(),
        panel.grid       = element_blank(),
        strip.background = element_rect(fill= NA),
        strip.placement = "outside",
        legend.background = element_rect(fill=NA,
                                         size=0.5),
        plot.margin = unit(c(0.5,1,0.5,1), "cm"),
        legend.position = "bottom")



dev.off()

###################################BD1######
new_BD1 = bind_cols(taiji_res[,30],taiji_res[,32],taiji_res[,34],taiji_res[,36],taiji_res[,38],taiji_res[,40],
                    taiji_res[,54],taiji_res[,56],taiji_res[,58],taiji_res[,60],taiji_res[,62],taiji_res[,64])
colnames(new_BD1) = c('GroupB_1','GroupB_2','GroupB_3','GroupB_4','GroupB_5','GroupB_6',
                      'GroupD_1','GroupD_2','GroupD_3','GroupD_4','GroupD_5','GroupD_6')
new_BD1$TFs =taiji_res$X
zscore_BD1 = new_BD1
for (i in 1:(ncol(new_BD1)-1)) {
  zscore_BD1[,i] = scale(new_BD1[,i])
}
zscore_BD1$mean = rowMeans(zscore_BD1[,1:12])
zscore_BD1$sd = rowSds(as.matrix(zscore_BD1[,1:12]))
tmp = data.frame((as.matrix(zscore_BD1[,1:12])-zscore_BD1$mean)/zscore_BD1$sd)
tmp$TFs = taiji_res$X
# setwd("~/Documents/Virginia Tech/Project_3_MIA mouse/taiji_new/ABCD/BD1")
# TFs_B = read.csv('potential_key_TFs_BD1_GroupB1.csv') 
# TFs_D = read.csv('potential_key_TFs_BD1_GroupD1.csv')
# TFs = union(TFs_B$taiji_res.X,TFs_D$taiji_res.X)
# setwd("~/Documents/Virginia Tech/Project_3_MIA mouse/taiji/ABCD/BD")
diff_TF = read.csv('Diff_TFs_BD1.csv')
# diff_TF = diff_TF[which(diff_TF$TFs %in% TFs),]
increase_TFs = diff_TF[which(diff_TF$diff > 0),]
decreade_TFs =diff_TF[which(diff_TF$diff < 0),]
# up_TF = TFs[which(TFs$V1 %in% increase_TFs$TFs),]
# down_TF = TFs[which(TFs$V1 %in% decreade_TFs$TFs),]
# write.table(up_TF,'up_TFs.txt',row.names = FALSE,quote = FALSE)
# write.table(down_TF,'down_TF.txt',row.names = FALSE,quote = FALSE)
TFs_matrix_up = tmp[which(tmp$TFs %in% increase_TFs$TFs),]
TFs_matrix_down = tmp[which(tmp$TFs %in% decreade_TFs$TFs),]
TFs_matrix = bind_rows(TFs_matrix_up,TFs_matrix_down)
new = as.matrix(TFs_matrix[,1:12])
rownames(new)= TFs_matrix$TFs
melted_BD1 = melt(new)

melted_BD1$Group <- NA
melted_BD1$Group[1:66] <- "GroupB"
melted_BD1$Group[67:132] <- "GroupD1"
pdf(file = "heatmap_BD1_key.pdf",width = 5,height = 3)
ggplot(data = melted_BD1, aes(x=Var2, y=Var1, fill= value)) + 
  geom_tile() + scale_fill_viridis(name = "PRP z-score") +
  facet_grid2(.~Group, scales = "free_x", space = "free_x", switch = "x")+
  guides(fill = guide_colourbar(title = "PRP z-score"), color= "black")+
  theme(axis.title.x     = element_blank(),
        axis.title.y     = element_blank(),
        axis.text.x      = element_blank(),
        panel.spacing.x  = unit(0, "mm"),
        panel.spacing.y  = unit(0, "mm"),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background = element_blank(),
        panel.grid       = element_blank(),
        strip.background = element_rect(fill= NA),
        strip.placement = "outside",
        legend.background = element_rect(fill=NA,
                                         size=0.5),
        plot.margin = unit(c(0.5,1,0.5,1), "cm"),
        legend.position = "bottom")



dev.off()
###################################BD2######
new_BD2 = bind_cols(taiji_res[,30],taiji_res[,32],taiji_res[,34],taiji_res[,36],taiji_res[,38],taiji_res[,40],
                    taiji_res[,66],taiji_res[,68],taiji_res[,70])
colnames(new_BD2) = c('GroupB_1','GroupB_2','GroupB_3','GroupB_4','GroupB_5','GroupB_6',
                      'GroupD_7','GroupD_8','GroupD_9')
new_BD2$TFs =taiji_res$X
zscore_BD2 = new_BD2
for (i in 1:(ncol(new_BD2)-1)) {
  zscore_BD2[,i] = scale(new_BD2[,i])
}
zscore_BD2$mean = rowMeans(zscore_BD2[,1:9])
zscore_BD2$sd = rowSds(as.matrix(zscore_BD2[,1:9]))
tmp = data.frame((as.matrix(zscore_BD2[,1:9])-zscore_BD2$mean)/zscore_BD2$sd)
tmp$TFs = taiji_res$X
# setwd("~/Documents/Virginia Tech/Project_3_MIA mouse/taiji_new/ABCD/BD2")
# TFs_B = read.csv('potential_key_TFs_BD2_GroupB2.csv') 
# TFs_D = read.csv('potential_key_TFs_BD2_GroupD2.csv')
# TFs = union(TFs_B$taiji_res.X,TFs_D$taiji_res.X)
# setwd("~/Documents/Virginia Tech/Project_3_MIA mouse/taiji/ABCD/BD/")
diff_TF = read.csv('Diff_TFs_BD2.csv')
# diff_TF = diff_TF[which(diff_TF$TFs %in% TFs),]
increase_TFs = diff_TF[which(diff_TF$diff > 0),]
decreade_TFs =diff_TF[which(diff_TF$diff < 0),]
# up_TF = TFs[which(TFs$V1 %in% increase_TFs$TFs),]
# down_TF = TFs[which(TFs$V1 %in% decreade_TFs$TFs),]
# write.table(up_TF,'up_TFs.txt',row.names = FALSE,quote = FALSE)
# write.table(down_TF,'down_TF.txt',row.names = FALSE,quote = FALSE)
TFs_matrix_up = tmp[which(tmp$TFs %in% increase_TFs$TFs),]
TFs_matrix_down = tmp[which(tmp$TFs %in% decreade_TFs$TFs),]
TFs_matrix = bind_rows(TFs_matrix_up,TFs_matrix_down)
new = as.matrix(TFs_matrix[,1:9])
rownames(new)= TFs_matrix$TFs
melted_BD2 = melt(new)

melted_BD2$Group <- NA
melted_BD2$Group[1:42] <- "GroupB"
melted_BD2$Group[43:63] <- "GroupD2"
pdf(file = "heatmap_BD2_key.pdf",width = 5,height = 3)
ggplot(data = melted_BD2, aes(x=Var2, y=Var1, fill= value)) + 
  geom_tile() + scale_fill_viridis(name = "PRP z-score") +
  facet_grid2(.~Group, scales = "free_x", space = "free_x", switch = "x")+
  guides(fill = guide_colourbar(title = "PRP z-score"), color= "black")+
  theme(axis.title.x     = element_blank(),
        axis.title.y     = element_blank(),
        axis.text.x      = element_blank(),
        panel.spacing.x  = unit(0, "mm"),
        panel.spacing.y  = unit(0, "mm"),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background = element_blank(),
        panel.grid       = element_blank(),
        strip.background = element_rect(fill= NA),
        strip.placement = "outside",
        legend.background = element_rect(fill=NA,
                                         size=0.5),
        plot.margin = unit(c(0.5,1,0.5,1), "cm"),
        legend.position = "bottom")



dev.off()
###################################BD3######
new_BD3 = bind_cols(taiji_res[,30],taiji_res[,32],taiji_res[,34],taiji_res[,36],taiji_res[,38],taiji_res[,40],
                    taiji_res[,72],taiji_res[,74],taiji_res[,76],taiji_res[,77],taiji_res[,79])
colnames(new_BD3) = c('GroupB_1','GroupB_2','GroupB_3','GroupB_4','GroupB_5','GroupB_6',
                      'GroupD_11','GroupD_12','GroupD_13','GroupD_14','GroupD_15')
new_BD3$TFs =taiji_res$X
zscore_BD3 = new_BD3
for (i in 1:(ncol(new_BD3)-1)) {
  zscore_BD3[,i] = scale(new_BD3[,i])
}
zscore_BD3$mean = rowMeans(zscore_BD3[,1:11])
zscore_BD3$sd = rowSds(as.matrix(zscore_BD3[,1:11]))
tmp = data.frame((as.matrix(zscore_BD3[,1:11])-zscore_BD3$mean)/zscore_BD3$sd)
tmp$TFs = taiji_res$X
# setwd("~/Documents/Virginia Tech/Project_3_MIA mouse/taiji_new/ABCD/BD3")
# TFs_B = read.csv('potential_key_TFs_BD3_GroupB3.csv') 
# TFs_D = read.csv('potential_key_TFs_BD3_GroupD3.csv')
# TFs = union(TFs_B$taiji_res.X,TFs_D$taiji_res.X)
# setwd("~/Documents/Virginia Tech/Project_3_MIA mouse/taiji/ABCD/BD/")
diff_TF = read.csv('Diff_TFs_BD3.csv')
# diff_TF = diff_TF[which(diff_TF$TFs %in% TFs),]
increase_TFs = diff_TF[which(diff_TF$diff > 0),]
decreade_TFs =diff_TF[which(diff_TF$diff < 0),]
# up_TF = TFs[which(TFs$V1 %in% increase_TFs$TFs),]
# down_TF = TFs[which(TFs$V1 %in% decreade_TFs$TFs),]
# write.table(up_TF,'up_TFs.txt',row.names = FALSE,quote = FALSE)
# write.table(down_TF,'down_TF.txt',row.names = FALSE,quote = FALSE)
TFs_matrix_up = tmp[which(tmp$TFs %in% increase_TFs$TFs),]
TFs_matrix_down = tmp[which(tmp$TFs %in% decreade_TFs$TFs),]
TFs_matrix = bind_rows(TFs_matrix_up,TFs_matrix_down)
new = as.matrix(TFs_matrix[,1:11])
rownames(new)= TFs_matrix$TFs
melted_BD3 = melt(new)

melted_BD3$Group <- NA
melted_BD3$Group[1:48] <- "GroupB"
melted_BD3$Group[49:88] <- "GroupD3"
pdf(file = "heatmap_BD3_key.pdf",width = 5,height = 3)
ggplot(data = melted_BD3, aes(x=Var2, y=Var1, fill= value)) +
  geom_tile() + scale_fill_viridis(name = "PRP z-score") +
  facet_grid2(.~Group, scales = "free_x", space = "free_x", switch = "x")+
  guides(fill = guide_colourbar(title = "PRP z-score"), color= "black")+
  theme(axis.title.x     = element_blank(),
        axis.title.y     = element_blank(),
        axis.text.x      = element_blank(),
        panel.spacing.x  = unit(0, "mm"),
        panel.spacing.y  = unit(0, "mm"),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background = element_blank(),
        panel.grid       = element_blank(),
        strip.background = element_rect(fill= NA),
        strip.placement = "outside",
        legend.background = element_rect(fill=NA,
                                         size=0.5),
        plot.margin = unit(c(0.5,1,0.5,1), "cm"),
        legend.position = "bottom")



dev.off()

###################################A1E######
new_A1E = bind_cols(taiji_res[,2],taiji_res[,4],taiji_res[,6],taiji_res[,7],taiji_res[,9],taiji_res[,11],
                    taiji_res_E[,2],taiji_res_E[,4],taiji_res_E[,6],taiji_res_E[,8])
colnames(new_A1E) = c('GroupA_1','GroupA_2','GroupA_3','GroupA_4','GroupA_5','GroupA_6',
                      'GroupE_2','GroupE_4','GroupE_5','GroupE_6')
new_A1E$TFs =taiji_res$X
zscore_A1E = new_A1E
for (i in 1:(ncol(new_A1E)-1)) {
  zscore_A1E[,i] = scale(new_A1E[,i])
}
zscore_A1E$mean = rowMeans(zscore_A1E[,1:10])
zscore_A1E$sd = rowSds(as.matrix(zscore_A1E[,1:10]))
tmp = data.frame((as.matrix(zscore_A1E[,1:10])-zscore_A1E$mean)/zscore_A1E$sd)
tmp$TFs = taiji_res$X
# setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/taiji/output_from_taiji/")
# TFs_B = read.csv('potential_key_TFs_A1E_GroupB1.csv') 
# TFs_A = read.csv('potential_key_TFs_A1E_GroupA1.csv')
# TFs = union(TFs_B$taiji_res.X,TFs_A$taiji_res.X)
# setwd("~/Documents/Virginia Tech/Project_3_MIA mouse/taiji/ABCD/BA")
diff_TF = read.csv('Diff_TFs_A1E.csv')
# diff_TF = diff_TF[which(diff_TF$TFs %in% TFs),]
increase_TFs = diff_TF[which(diff_TF$diff > 0),]
decreade_TFs =diff_TF[which(diff_TF$diff < 0),]
# up_TF = TFs[which(TFs$V1 %in% increase_TFs$TFs),]
# down_TF = TFs[which(TFs$V1 %in% decreade_TFs$TFs),]
# write.table(up_TF,'up_TFs.txt',row.names = FALSE,quote = FALSE)
# write.table(down_TF,'down_TF.txt',row.names = FALSE,quote = FALSE)
TFs_matrix_up = tmp[which(tmp$TFs %in% increase_TFs$TFs),]
TFs_matrix_down = tmp[which(tmp$TFs %in% decreade_TFs$TFs),]
TFs_matrix = bind_rows(TFs_matrix_up,TFs_matrix_down)
new = as.matrix(TFs_matrix[,1:10])
rownames(new)= TFs_matrix$TFs
melted_A1E = melt(new)

melted_A1E$Group <- NA
melted_A1E$Group[1:48] <- "GroupA1"
melted_A1E$Group[49:80] <- "GroupE"
pdf(file = "heatmap_A1E_key.pdf",width = 5,height = 3)
ggplot(data = melted_A1E, aes(x=Var2, y=Var1, fill= value)) +
  geom_tile() + scale_fill_viridis(name = "PRP z-score") +
  facet_grid2(.~Group, scales = "free_x", space = "free_x", switch = "x")+
  guides(fill = guide_colourbar(title = "PRP z-score"), color= "black")+
  theme(axis.title.x     = element_blank(),
        axis.title.y     = element_blank(),
        axis.text.x      = element_blank(),
        panel.spacing.x  = unit(0, "mm"),
        panel.spacing.y  = unit(0, "mm"),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background = element_blank(),
        panel.grid       = element_blank(),
        strip.background = element_rect(fill= NA),
        strip.placement = "outside",
        legend.background = element_rect(fill=NA,
                                         size=0.5),
        legend.key = element_rect(colour = "black", size=0.5),
        plot.margin = unit(c(0.5,1,0.5,1), "cm"),
        legend.position = "bottom")

dev.off()

###################################A2E######
new_A2E = bind_cols(taiji_res[,13],taiji_res[,15],taiji_res[,17], taiji_res_E[,2],taiji_res_E[,4],taiji_res_E[,6],taiji_res_E[,8])
colnames(new_A2E) = c('GroupA_7','GroupA_8','GroupA_9','GroupE_2','GroupE_4','GroupE_5','GroupE_6')
new_A2E$TFs =taiji_res$X
zscore_A2E = new_A2E
for (i in 1:(ncol(new_A2E)-1)) {
  zscore_A2E[,i] = scale(new_A2E[,i])
}
zscore_A2E$mean = rowMeans(zscore_A2E[,1:7])
zscore_A2E$sd = rowSds(as.matrix(zscore_A2E[,1:7]))
tmp = data.frame((as.matrix(zscore_A2E[,1:7])-zscore_A2E$mean)/zscore_A2E$sd)
tmp$TFs = taiji_res$X
# setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/taiji/output_from_taiji/")
# TFs_B = read.csv('potential_key_TFs_A2E_GroupB1.csv') 
# TFs_A = read.csv('potential_key_TFs_A2E_GroupA1.csv')
# TFs = union(TFs_B$taiji_res.X,TFs_A$taiji_res.X)
# setwd("~/Documents/Virginia Tech/Project_3_MIA mouse/taiji/ABCD/BA")
diff_TF = read.csv('Diff_TFs_A2E.csv')
# diff_TF = diff_TF[which(diff_TF$TFs %in% TFs),]
increase_TFs = diff_TF[which(diff_TF$diff > 0),]
decreade_TFs =diff_TF[which(diff_TF$diff < 0),]
# up_TF = TFs[which(TFs$V1 %in% increase_TFs$TFs),]
# down_TF = TFs[which(TFs$V1 %in% decreade_TFs$TFs),]
# write.table(up_TF,'up_TFs.txt',row.names = FALSE,quote = FALSE)
# write.table(down_TF,'down_TF.txt',row.names = FALSE,quote = FALSE)
TFs_matrix_up = tmp[which(tmp$TFs %in% increase_TFs$TFs),]
TFs_matrix_down = tmp[which(tmp$TFs %in% decreade_TFs$TFs),]
TFs_matrix = bind_rows(TFs_matrix_up,TFs_matrix_down)
new = as.matrix(TFs_matrix[,1:7])
rownames(new)= TFs_matrix$TFs
melted_A2E = melt(new)

melted_A2E$Group <- NA
melted_A2E$Group[1:24] <- "GroupA2"
melted_A2E$Group[25:56] <- "GroupE"
pdf(file = "heatmap_A2E_key.pdf",width = 5,height = 3)
ggplot(data = melted_A2E, aes(x=Var2, y=Var1, fill= value)) +
  geom_tile() + scale_fill_viridis(name = "PRP z-score") +
  facet_grid2(.~Group, scales = "free_x", space = "free_x", switch = "x")+
  guides(fill = guide_colourbar(title = "PRP z-score"), color= "black")+
  theme(axis.title.x     = element_blank(),
        axis.title.y     = element_blank(),
        axis.text.x      = element_blank(),
        panel.spacing.x  = unit(0, "mm"),
        panel.spacing.y  = unit(0, "mm"),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background = element_blank(),
        panel.grid       = element_blank(),
        strip.background = element_rect(fill= NA),
        strip.placement = "outside",
        legend.background = element_rect(fill=NA,
                                         size=0.5),
        legend.key = element_rect(colour = "black", size=0.5),
        plot.margin = unit(c(0.5,1,0.5,1), "cm"),
        legend.position = "bottom")

dev.off()
###################################A3E######
new_A3E = bind_cols(taiji_res[,19],taiji_res[,20],taiji_res[,22],taiji_res[,24],
                    taiji_res[,26],taiji_res[,28], taiji_res_E[,2],taiji_res_E[,4],taiji_res_E[,6],taiji_res_E[,8])
colnames(new_A3E) = c('GroupA_10','GroupA_11','GroupA_12','GroupA_13',
                      'GroupA_14','GroupA_15','GroupE_2','GroupE_4','GroupE_5','GroupE_6')
new_A3E$TFs =taiji_res$X
zscore_A3E = new_A3E
for (i in 1:(ncol(new_A3E)-1)) {
  zscore_A3E[,i] = scale(new_A3E[,i])
}
zscore_A3E$mean = rowMeans(zscore_A3E[,1:10])
zscore_A3E$sd = rowSds(as.matrix(zscore_A3E[,1:10]))
tmp = data.frame((as.matrix(zscore_A3E[,1:10])-zscore_A3E$mean)/zscore_A3E$sd)
tmp$TFs = taiji_res$X
# setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/taiji/output_from_taiji/")
# TFs_B = read.csv('potential_key_TFs_A3E_GroupB1.csv') 
# TFs_A = read.csv('potential_key_TFs_A3E_GroupA1.csv')
# TFs = union(TFs_B$taiji_res.X,TFs_A$taiji_res.X)
# setwd("~/Documents/Virginia Tech/Project_3_MIA mouse/taiji/ABCD/BA")
diff_TF = read.csv('Diff_TFs_A3E.csv')
# diff_TF = diff_TF[which(diff_TF$TFs %in% TFs),]
increase_TFs = diff_TF[which(diff_TF$diff > 0),]
decreade_TFs =diff_TF[which(diff_TF$diff < 0),]
# up_TF = TFs[which(TFs$V1 %in% increase_TFs$TFs),]
# down_TF = TFs[which(TFs$V1 %in% decreade_TFs$TFs),]
# write.table(up_TF,'up_TFs.txt',row.names = FALSE,quote = FALSE)
# write.table(down_TF,'down_TF.txt',row.names = FALSE,quote = FALSE)
TFs_matrix_up = tmp[which(tmp$TFs %in% increase_TFs$TFs),]
TFs_matrix_down = tmp[which(tmp$TFs %in% decreade_TFs$TFs),]
TFs_matrix = bind_rows(TFs_matrix_up,TFs_matrix_down)
new = as.matrix(TFs_matrix[,1:10])
rownames(new)= TFs_matrix$TFs
melted_A3E = melt(new)

melted_A3E$Group <- NA
melted_A3E$Group[1:42] <- "GroupA3"
melted_A3E$Group[43:70] <- "GroupE"
pdf(file = "heatmap_A3E_key.pdf",width = 5,height = 3)
ggplot(data = melted_A3E, aes(x=Var2, y=Var1, fill= value)) +
  geom_tile() + scale_fill_viridis(name = "PRP z-score") +
  facet_grid2(.~Group, scales = "free_x", space = "free_x", switch = "x")+
  guides(fill = guide_colourbar(title = "PRP z-score"), color= "black")+
  theme(axis.title.x     = element_blank(),
        axis.title.y     = element_blank(),
        axis.text.x      = element_blank(),
        panel.spacing.x  = unit(0, "mm"),
        panel.spacing.y  = unit(0, "mm"),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background = element_blank(),
        panel.grid       = element_blank(),
        strip.background = element_rect(fill= NA),
        strip.placement = "outside",
        legend.background = element_rect(fill=NA,
                                         size=0.5),
        legend.key = element_rect(colour = "black", size=0.5),
        plot.margin = unit(c(0.5,1,0.5,1), "cm"),
        legend.position = "bottom")

dev.off()

