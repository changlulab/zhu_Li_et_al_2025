###########prep#############
rm(list = ls())
options(stringsAsFactors = F)
library(dplyr)
library(MASS)
library(fitdistrplus)
library(reshape2)
library(ggplot2)
library(DESeq2)
library(matrixStats)
library(ggpubr)
library(ggh4x)
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/taiji/output_from_taiji")
taiji_res = read.csv('GeneRanks_EFG.tsv',sep = '\t',header = T)
###################################FG######
new_FG = bind_cols(taiji_res[,10],taiji_res[,12],taiji_res[,14],taiji_res[,16],taiji_res[,18],
                   taiji_res[,20],taiji_res[,22],taiji_res[,24],taiji_res[,26],taiji_res[,28], taiji_res[,30])
colnames(new_FG) = c('GroupF_1','GroupF_2','GroupF_3','GroupF_5','GroupF_6',
                     'GroupG_1','GroupG_2','GroupG_3','GroupG_4','GroupG_5','GroupG_6')
new_FG$TFs =taiji_res$X
zscore_FG = new_FG
for (i in 1:(ncol(new_FG)-1)) {
  zscore_FG[,i] = scale(new_FG[,i])
}
zscore_FG$mean = rowMeans(zscore_FG[,1:11])
zscore_FG$sd = rowSds(as.matrix(zscore_FG[,1:11]))
tmp = data.frame((as.matrix(zscore_FG[,1:11])-zscore_FG$mean)/zscore_FG$sd)
tmp$TFs = taiji_res$X
# setwd("~/Documents/Virginia Tech/Project_3_MIA mouse/taiji_new/EFG/FG")
# TFs_F = read.csv('potential_key_TFs_GroupF.csv') 
# TFs_G = read.csv('potential_key_TFs_GroupG.csv')
# TFs = union(TFs_F$taiji_res.X,TFs_G$taiji_res.X)
# setwd("~/Documents/Virginia Tech/Project_3_MIA mouse/taiji/EFG/FG/")
diff_TF = read.csv('Diff_TFs_FG.csv')
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
melted_FG = melt(new)

melted_FG$Group <- NA
melted_FG$Group[1:75] <- "GroupF"
melted_FG$Group[76:165] <- "GroupG"
pdf(file = "heatmap_FG_key.pdf",width = 5,height = 4)
ggplot(data = melted_FG, aes(x=Var2, y=Var1, fill= value)) +
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
###################################FE######
new_FE = bind_cols(taiji_res[,10],taiji_res[,12],taiji_res[,14],taiji_res[,16],taiji_res[,18],
                   taiji_res[,2],taiji_res[,4],taiji_res[,6],taiji_res[,8])
colnames(new_FE) = c('GroupF_1','GroupF_2','GroupF_3','GroupF_4','GroupF_5',
                     'GroupE_2','GroupE_3','GroupE_4','GroupE_6')
new_FE$TFs =taiji_res$X
zscore_FE = new_FE
for (i in 1:(ncol(new_FE)-1)) {
  zscore_FE[,i] = scale(new_FE[,i])
}
zscore_FE$mean = rowMeans(zscore_FE[,1:9])
zscore_FE$sd = rowSds(as.matrix(zscore_FE[,1:9]))
tmp = data.frame((as.matrix(zscore_FE[,1:9])-zscore_FE$mean)/zscore_FE$sd)
tmp$TFs = taiji_res$X
# setwd("~/Documents/Virginia Tech/Project_3_MIA mouse/taiji_new/EFG/FE")
# TFs_F = read.csv('potential_key_TFs_GroupF.csv') 
# TFs_E = read.csv('potential_key_TFs_GroupE.csv')
# TFs = union(TFs_F$taiji_res.X,TFs_E$taiji_res.X)
# setwd("~/Documents/Virginia Tech/Project_3_MIA mouse/taiji/EFG/FE/")
diff_TF = read.csv('Diff_TFs_FE.csv')
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
melted_FE = melt(new)

melted_FE$Group <- NA
melted_FE$Group[1:60] <- "GroupF"
melted_FE$Group[61:108] <- "GroupE"
pdf(file = "heatmap_FE_key.pdf",width = 5,height = 4)
ggplot(data = melted_FE, aes(x=Var2, y=Var1, fill= value)) +
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

###################################EG######
new_EG = bind_cols(taiji_res[,2],taiji_res[,4],taiji_res[,6],taiji_res[,8],
                   taiji_res[,20],taiji_res[,22],taiji_res[,24],taiji_res[,26],taiji_res[,28])
colnames(new_EG) = c('GroupE_2','GroupE_4','GroupE_5','GroupE_6',
                     'GroupG_1','GroupG_2','GroupG_4','GroupG_5','GroupG_6')
new_EG$TFs =taiji_res$X
zscore_EG = new_EG
for (i in 1:(ncol(new_EG)-1)) {
  zscore_EG[,i] = scale(new_EG[,i])
}
zscore_EG$mean = rowMeans(zscore_EG[,1:9])
zscore_EG$sd = rowSds(as.matrix(zscore_EG[,1:9]))
tmp = data.frame((as.matrix(zscore_EG[,1:9])-zscore_EG$mean)/zscore_EG$sd)
tmp$TFs = taiji_res$X
# setwd("~/Documents/Virginia Tech/Project_3_MIA mouse/taiji_new/EFG/EG")
# TFs_E = read.csv('potential_key_TFs_GroupE.csv') 
# TFs_G = read.csv('potential_key_TFs_GroupG.csv')
# TFs = union(TFs_E$taiji_res.X,TFs_G$taiji_res.X)
# setwd("~/Documents/Virginia Tech/Project_3_MIA mouse/taiji/EFG/EG/")
diff_TF = read.csv('Diff_TFs_EG.csv')
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
melted_EG = melt(new)

melted_EG$Group <- NA
melted_EG$Group[1:72] <- "GroupE"
melted_EG$Group[73:162] <- "GroupG"
pdf(file = "heatmap_EG_key.pdf",width = 5,height = 4)
ggplot(data = melted_EG, aes(x=Var2, y=Var1, fill= value)) +
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
