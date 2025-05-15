##########preparation#########
rm(list = ls())
options(stringsAsFactors = F)
library(dplyr)
library(MASS)
library(fitdistrplus)
library(reshape2)
library(ggplot2)
library(DESeq2)
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/taiji/output_from_taiji")
# taiji_res = read.csv('GeneRanks_EFG.tsv',sep = '\t',header = T)
taiji_res = read.csv('GeneRanks_EFG_new.tsv',sep = '\t',header = T)
##################FG#####################
new_FG = bind_cols(taiji_res[,10],taiji_res[,12],taiji_res[,14],taiji_res[,16],taiji_res[,18],
                   taiji_res[,20],taiji_res[,22],taiji_res[,24],taiji_res[,26],taiji_res[,28], taiji_res[,30])
colnames(new_FG) = c('GroupF_1','GroupF_2','GroupF_3','GroupF_5','GroupF_6',
                     'GroupG_1','GroupG_2','GroupG_3','GroupG_4','GroupG_5','GroupG_6')
new_FG$TFs =taiji_res$X
row.names(new_FG) = taiji_res$X
aver_GroupF = rowMeans(new_FG[,1:5])
aver_GroupG = rowMeans(new_FG[,6:11])
#####################potential key TFs#######################
fit_GroupF = fitdistr(aver_GroupF,'normal')
para_GroupF = fit_GroupF$estimate
pvalue_GroupF = pnorm(aver_GroupF,mean = para_GroupF[1], sd = para_GroupF[2], lower.tail = F)
FDR_GroupF = p.adjust(pvalue_GroupF,method = "fdr", n = length(pvalue_GroupF))
p_GroupF = data.frame(aver_GroupF,pvalue_GroupF,FDR_GroupF,taiji_res$X)
key_GroupF = p_GroupF[which(p_GroupF$FDR_GroupF < 0.05),]
write.csv(key_GroupF,'potential_key_TFs_GroupF.csv',row.names = FALSE)

fit_GroupG = fitdistr(aver_GroupG,'normal')
para_GroupG = fit_GroupG$estimate
pvalue_GroupG = pnorm(aver_GroupG,mean = para_GroupG[1], sd = para_GroupG[2], lower.tail = F)
FDR_GroupG = p.adjust(pvalue_GroupG,method = "fdr", n = length(pvalue_GroupG))
p_GroupG = data.frame(aver_GroupG,pvalue_GroupG,FDR_GroupG,taiji_res$X)
key_GroupG = p_GroupG[which(p_GroupG$FDR_GroupG < 0.05),]
write.csv(key_GroupG,'potential_key_TFs_GroupG.csv',row.names = FALSE)
#############################differential#######################
tmp = scale(data.frame(aver_GroupF,aver_GroupG))
nor_mat_FG = data.frame(tmp)
row.names(nor_mat_FG) = taiji_res$X
colnames(nor_mat_FG) = c('GroupF','GroupG')
nor_mat_FG$diff = nor_mat_FG$GroupF - nor_mat_FG$GroupG
fit_diff_FG = fitdistr(abs(nor_mat_FG$diff),'normal')
para_diff_FG = fit_diff_FG$estimate
pvalue_diff_FG = pnorm(abs(nor_mat_FG$diff), mean = para_diff_FG[1], sd = para_diff_FG[2],lower.tail = F)
FDR_diff_FG = p.adjust(pvalue_diff_FG,method = "fdr", n = length(pvalue_diff_FG))
nor_mat_FG$pvalue = pvalue_diff_FG
nor_mat_FG$FDR = FDR_diff_FG
nor_mat_FG$TFs = row.names(nor_mat_FG)
diff_TFs_FG = nor_mat_FG[which(nor_mat_FG$FDR < 0.05),]
new_matrix_FG = new_FG[which(new_FG$TFs %in% diff_TFs_FG$TFs),]
final_diff_TFs_FG = inner_join(diff_TFs_FG,new_matrix_FG,by = 'TFs')
write.csv(final_diff_TFs_FG,'Diff_TFs_FG.csv',row.names = FALSE)
#################FE############
new_FE = bind_cols(taiji_res[,10],taiji_res[,12],taiji_res[,14],taiji_res[,16],taiji_res[,18],
                    taiji_res[,2],taiji_res[,4],taiji_res[,6],taiji_res[,8])
colnames(new_FE) = c('GroupF_1','GroupF_2','GroupF_3','GroupF_5','GroupF_6',
                      'GroupE_2','GroupE_4','GroupE_5','GroupE_6')
new_FE$TFs =taiji_res$X
row.names(new_FE) = taiji_res$X
aver_GroupF = rowMeans(new_FE[,1:5])
aver_GroupE = rowMeans(new_FE[,6:9])
#####################potential key TFs#######################
fit_GroupF = fitdistr(aver_GroupF,'normal')
para_GroupF = fit_GroupF$estimate
pvalue_GroupF = pnorm(aver_GroupF,mean = para_GroupF[1], sd = para_GroupF[2], lower.tail = F)
FDR_GroupF = p.adjust(pvalue_GroupF,method = "fdr", n = length(pvalue_GroupF))
p_GroupF = data.frame(aver_GroupF,pvalue_GroupF,FDR_GroupF,taiji_res$X)
key_GroupF = p_GroupF[which(p_GroupF$FDR_GroupF < 0.05),]
write.csv(key_GroupF,'potential_key_TFs_FE_GroupF.csv',row.names = FALSE)

fit_GroupE = fitdistr(aver_GroupE,'normal')
para_GroupE = fit_GroupE$estimate
pvalue_GroupE = pnorm(aver_GroupE,mean = para_GroupE[1], sd = para_GroupE[2], lower.tail = F)
FDR_GroupE = p.adjust(pvalue_GroupE,method = "fdr", n = length(pvalue_GroupE))
p_GroupE = data.frame(aver_GroupE,pvalue_GroupE,FDR_GroupE,taiji_res$X)
key_GroupE = p_GroupE[which(p_GroupE$FDR_GroupE < 0.05),]
write.csv(key_GroupE,'potential_key_TFs_FE_GroupE.csv',row.names = FALSE)
#############################differential#######################
tmp = scale(data.frame(aver_GroupF,aver_GroupE))
nor_mat_FE = data.frame(tmp)
row.names(nor_mat_FE) = taiji_res$X
colnames(nor_mat_FE) = c('GroupF','GrpupE')
nor_mat_FE$diff = nor_mat_FE$GroupF - nor_mat_FE$GrpupE
fit_diff_FE = fitdistr(abs(nor_mat_FE$diff),'normal')
para_diff_FE = fit_diff_FE$estimate
pvalue_diff_FE = pnorm(abs(nor_mat_FE$diff), mean = para_diff_FE[1], sd = para_diff_FE[2],lower.tail = F)
FDR_diff_FE = p.adjust(pvalue_diff_FE,method = "fdr", n = length(pvalue_diff_FE))
nor_mat_FE$pvalue = pvalue_diff_FE
nor_mat_FE$FDR = FDR_diff_FE
nor_mat_FE$TFs = row.names(nor_mat_FE)
diff_TFs_FE = nor_mat_FE[which(nor_mat_FE$FDR < 0.05),]
new_matrix_FE = new_FE[which(new_FE$TFs %in% diff_TFs_FE$TFs),]
final_diff_TFs_FE = inner_join(diff_TFs_FE,new_matrix_FE,by = 'TFs')
write.csv(final_diff_TFs_FE,'Diff_TFs_FE.csv',row.names = FALSE)
#################EG############
new_EG = bind_cols(taiji_res[,2],taiji_res[,4],taiji_res[,6],taiji_res[,8],
                   taiji_res[,20],taiji_res[,22],taiji_res[,24],taiji_res[,26],taiji_res[,28], taiji_res[,30])
colnames(new_EG) = c('GroupE_2','GroupE_4','GroupE_5','GroupE_6',
                     'GroupG_1','GroupG_2','GroupG_3','GroupG_4','GroupG_5','GroupG_6')
new_EG$TFs = taiji_res$X
row.names(new_EG) = taiji_res$X
aver_GroupE = rowMeans(new_EG[,1:4])
aver_GroupG = rowMeans(new_EG[,5:10])
#####################potential key TFs#######################
fit_GroupE = fitdistr(aver_GroupE,'normal')
para_GroupE = fit_GroupE$estimate
pvalue_GroupE = pnorm(aver_GroupE,mean = para_GroupE[1], sd = para_GroupE[2], lower.tail = F)
FDR_GroupE = p.adjust(pvalue_GroupE,method = "fdr", n = length(pvalue_GroupE))
p_GroupE = data.frame(aver_GroupE,pvalue_GroupE,FDR_GroupE,taiji_res$X)
key_GroupE = p_GroupE[which(p_GroupE$FDR_GroupE < 0.05),]
write.csv(key_GroupE,'potential_key_TFs_GroupE.csv',row.names = FALSE)

fit_GroupG = fitdistr(aver_GroupG,'normal')
para_GroupG = fit_GroupG$estimate
pvalue_GroupG = pnorm(aver_GroupG,mean = para_GroupG[1], sd = para_GroupG[2], lower.tail = F)
FDR_GroupG = p.adjust(pvalue_GroupG,method = "fdr", n = length(pvalue_GroupG))
p_GroupG = data.frame(aver_GroupG,pvalue_GroupG,FDR_GroupG,taiji_res$X)
key_GroupG = p_GroupG[which(p_GroupG$FDR_GroupG < 0.05),]
write.csv(key_GroupG,'potential_key_TFs_GroupG.csv',row.names = FALSE)
#############################differential#######################
tmp = scale(data.frame(aver_GroupE,aver_GroupG))
nor_mat_EG = data.frame(tmp)
row.names(nor_mat_EG) = taiji_res$X
colnames(nor_mat_EG) = c('GroupE','GroupG')
nor_mat_EG$diff = nor_mat_EG$GroupE - nor_mat_EG$GroupG
fit_diff_EG = fitdistr(abs(nor_mat_EG$diff),'normal')
para_diff_EG = fit_diff_EG$estimate
pvalue_diff_EG = pnorm(abs(nor_mat_EG$diff), mean = para_diff_EG[1], sd = para_diff_EG[2],lower.tail = F)
FDR_diff_EG = p.adjust(pvalue_diff_EG,method = "fdr", n = length(pvalue_diff_EG))
nor_mat_EG$pvalue = pvalue_diff_EG
nor_mat_EG$FDR = FDR_diff_EG
nor_mat_EG$TFs = row.names(nor_mat_EG)
diff_TFs_EG = nor_mat_EG[which(nor_mat_EG$FDR < 0.05),]
new_matrix_EG = new_EG[which(new_EG$TFs %in% diff_TFs_EG$TFs),]
final_diff_TFs_EG = inner_join(diff_TFs_EG,new_matrix_EG,by = 'TFs')
write.csv(final_diff_TFs_EG,'Diff_TFs_EG.csv',row.names = FALSE)
