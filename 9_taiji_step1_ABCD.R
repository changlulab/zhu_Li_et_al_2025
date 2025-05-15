rm(list = ls())
options(stringsAsFactors = F)
library(dplyr)
library(MASS)
library(fitdistrplus)
library(reshape2)
library(ggplot2)
library(DESeq2)
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/taiji/output_from_taiji/")
taiji_res = read.csv('GeneRanks_ABCD.tsv',sep = '\t',header = T)
taiji_res_E = read.csv('GeneRanks_EFG.tsv',sep = '\t',header = T)
##################BA1#####################
new_BA1 = bind_cols(taiji_res[,30],taiji_res[,32],taiji_res[,34],taiji_res[,36],taiji_res[,38],taiji_res[,40],
                    taiji_res[,2],taiji_res[,4],taiji_res[,6],taiji_res[,7],taiji_res[,9],taiji_res[,11])
colnames(new_BA1) = c('GroupB_1','GroupB_2','GroupB_3','GroupB_4','GroupB_5','GroupB_6',
                      'GroupA_1','GroupA_2','GroupA_3','GroupA_4','GroupA_5','GroupA_6')
new_BA1$TFs =taiji_res$X
row.names(new_BA1) = taiji_res$X
aver_GroupB1 = rowMeans(new_BA1[,1:6])
aver_GroupA1 = rowMeans(new_BA1[,7:12])
#####################potential key TFs#######################
fit_GroupB1 = fitdistr(aver_GroupB1,'normal')
para_GroupB1 = fit_GroupB1$estimate
pvalue_GroupB1 = pnorm(aver_GroupB1,mean = para_GroupB1[1], sd = para_GroupB1[2], lower.tail = F)
FDR_GroupB1 = p.adjust(pvalue_GroupB1,method = "fdr", n = length(pvalue_GroupB1))
p_GroupB1 = data.frame(aver_GroupB1,pvalue_GroupB1,FDR_GroupB1,taiji_res$X)
key_GroupB1 = p_GroupB1[which(p_GroupB1$FDR_GroupB1 < 0.05),]
write.csv(key_GroupB1,'potential_key_TFs_BA1_GroupB1.csv',row.names = FALSE)

fit_GroupA1 = fitdistr(aver_GroupA1,'normal')
para_GroupA1 = fit_GroupA1$estimate
pvalue_GroupA1 = pnorm(aver_GroupA1,mean = para_GroupA1[1], sd = para_GroupA1[2], lower.tail = F)
FDR_GroupA1 = p.adjust(pvalue_GroupA1,method = "fdr", n = length(pvalue_GroupA1))
p_GroupA1 = data.frame(aver_GroupA1,pvalue_GroupA1,FDR_GroupA1,taiji_res$X)
key_GroupA1 = p_GroupA1[which(p_GroupA1$FDR_GroupA1 < 0.05),]
write.csv(key_GroupA1,'potential_key_TFs_BA1_GroupA1.csv',row.names = FALSE)
#############################differential#######################
tmp = scale(data.frame(aver_GroupB1,aver_GroupA1))
nor_mat_BA1 = data.frame(tmp)
row.names(nor_mat_BA1) = taiji_res$X
colnames(nor_mat_BA1) = c('GroupB1','GroupA1')
nor_mat_BA1$diff = nor_mat_BA1$GroupB1 - nor_mat_BA1$GroupA1
fit_diff_BA1 = fitdistr(abs(nor_mat_BA1$diff),'normal')
para_diff_BA1 = fit_diff_BA1$estimate
pvalue_diff_BA1 = pnorm(abs(nor_mat_BA1$diff), mean = para_diff_BA1[1], sd = para_diff_BA1[2],lower.tail = F)
FDR_diff_BA1 = p.adjust(pvalue_diff_BA1,method = "fdr", n = length(pvalue_diff_BA1))
nor_mat_BA1$pvalue = pvalue_diff_BA1
nor_mat_BA1$FDR = FDR_diff_BA1
nor_mat_BA1$TFs = row.names(nor_mat_BA1)
diff_TFs_BA1 = nor_mat_BA1[which(nor_mat_BA1$FDR < 0.05),]
new_matrix_BA1 = new_BA1[which(new_BA1$TFs %in% diff_TFs_BA1$TFs),]
final_diff_TFs_BA1 = inner_join(diff_TFs_BA1,new_matrix_BA1,by = 'TFs')
write.csv(final_diff_TFs_BA1,'Diff_TFs_BA1.csv',row.names = FALSE)

##################BA2#####################
new_BA2 = bind_cols(taiji_res[,30],taiji_res[,32],taiji_res[,34],taiji_res[,36],taiji_res[,38],taiji_res[,40],
                    taiji_res[,13],taiji_res[,15],taiji_res[,17])
colnames(new_BA2) = c('GroupB_1','GroupB_2','GroupB_3','GroupB_4',
                      'GroupB_5','GroupB_6',
                      'GroupA_7','GroupA_8','GroupA_9')
new_BA2$TFs =taiji_res$X
row.names(new_BA2) = taiji_res$X
aver_GroupB2 = rowMeans(new_BA2[,1:6])
aver_GroupA2 = rowMeans(new_BA2[,7:9])
#####################potential key TFs#######################
fit_GroupB2 = fitdistr(aver_GroupB2,'normal')
para_GroupB2 = fit_GroupB2$estimate
pvalue_GroupB2 = pnorm(aver_GroupB2,mean = para_GroupB2[1], sd = para_GroupB2[2], lower.tail = F)
FDR_GroupB2 = p.adjust(pvalue_GroupB2,method = "fdr", n = length(pvalue_GroupB2))
p_GroupB2 = data.frame(aver_GroupB2,pvalue_GroupB2,FDR_GroupB2,taiji_res$X)
key_GroupB2 = p_GroupB2[which(p_GroupB2$FDR_GroupB2 < 0.05),]
write.csv(key_GroupB2,'potential_key_TFs_BA2_GroupB2.csv',row.names = FALSE)

fit_GroupA2 = fitdistr(aver_GroupA2,'normal')
para_GroupA2 = fit_GroupA2$estimate
pvalue_GroupA2 = pnorm(aver_GroupA2,mean = para_GroupA2[1], sd = para_GroupA2[2], lower.tail = F)
FDR_GroupA2 = p.adjust(pvalue_GroupA2,method = "fdr", n = length(pvalue_GroupA2))
p_GroupA2 = data.frame(aver_GroupA2,pvalue_GroupA2,FDR_GroupA2,taiji_res$X)
key_GroupA2 = p_GroupA2[which(p_GroupA2$FDR_GroupA2 < 0.05),]
write.csv(key_GroupA2,'potential_key_TFs_BA2_GroupA2.csv',row.names = FALSE)
#############################differential#######################
tmp = scale(data.frame(aver_GroupB2,aver_GroupA2))
nor_mat_BA2 = data.frame(tmp)
row.names(nor_mat_BA2) = taiji_res$X
colnames(nor_mat_BA2) = c('GroupB2','GroupA2')
nor_mat_BA2$diff = nor_mat_BA2$GroupB2 - nor_mat_BA2$GroupA2
fit_diff_BA2 = fitdistr(abs(nor_mat_BA2$diff),'normal')
para_diff_BA2 = fit_diff_BA2$estimate
pvalue_diff_BA2 = pnorm(abs(nor_mat_BA2$diff), mean = para_diff_BA2[1], sd = para_diff_BA2[2],lower.tail = F)
FDR_diff_BA2 = p.adjust(pvalue_diff_BA2,method = "fdr", n = length(pvalue_diff_BA2))
nor_mat_BA2$pvalue = pvalue_diff_BA2
nor_mat_BA2$FDR = FDR_diff_BA2
nor_mat_BA2$TFs = row.names(nor_mat_BA2)
diff_TFs_BA2 = nor_mat_BA2[which(nor_mat_BA2$FDR < 0.05),]
new_matrix_BA2 = new_BA2[which(new_BA2$TFs %in% diff_TFs_BA2$TFs),]
final_diff_TFs_BA2 = inner_join(diff_TFs_BA2,new_matrix_BA2,by = 'TFs')
write.csv(final_diff_TFs_BA2,'Diff_TFs_BA2.csv',row.names = FALSE)

##################BA3#####################
new_BA3 = bind_cols(taiji_res[,30],taiji_res[,32],taiji_res[,34],taiji_res[,36],taiji_res[,38],taiji_res[,40],
                    taiji_res[,19],taiji_res[,20],taiji_res[,22],taiji_res[,24],
                    taiji_res[,26],taiji_res[,28])
colnames(new_BA3) = c('GroupB_1','GroupB_2','GroupB_3','GroupB_4',
                      'GroupB_5','GroupB_6',
                      'GroupA_10','GroupA_11','GroupA_12','GroupA_13',
                      'GroupA_14','GroupA_15')
new_BA3$TFs =taiji_res$X
row.names(new_BA3) = taiji_res$X
aver_GroupB3 = rowMeans(new_BA3[,1:6])
aver_GroupA3 = rowMeans(new_BA3[,7:12])
#####################potential key TFs#######################
fit_GroupB3 = fitdistr(aver_GroupB3,'normal')
para_GroupB3 = fit_GroupB3$estimate
pvalue_GroupB3 = pnorm(aver_GroupB3,mean = para_GroupB3[1], sd = para_GroupB3[2], lower.tail = F)
FDR_GroupB3 = p.adjust(pvalue_GroupB3,method = "fdr", n = length(pvalue_GroupB3))
p_GroupB3 = data.frame(aver_GroupB3,pvalue_GroupB3,FDR_GroupB3,taiji_res$X)
key_GroupB3 = p_GroupB3[which(p_GroupB3$FDR_GroupB3 < 0.05),]
write.csv(key_GroupB3,'potential_key_TFs_BA3_GroupB3.csv',row.names = FALSE)

fit_GroupA3 = fitdistr(aver_GroupA3,'normal')
para_GroupA3 = fit_GroupA3$estimate
pvalue_GroupA3 = pnorm(aver_GroupA3,mean = para_GroupA3[1], sd = para_GroupA3[2], lower.tail = F)
FDR_GroupA3 = p.adjust(pvalue_GroupA3,method = "fdr", n = length(pvalue_GroupA3))
p_GroupA3 = data.frame(aver_GroupA3,pvalue_GroupA3,FDR_GroupA3,taiji_res$X)
key_GroupA3 = p_GroupA3[which(p_GroupA3$FDR_GroupA3 < 0.05),]
write.csv(key_GroupA3,'potential_key_TFs_BA3_GroupA3.csv',row.names = FALSE)
#############################differential#######################
tmp = scale(data.frame(aver_GroupB3,aver_GroupA3))
nor_mat_BA3 = data.frame(tmp)
row.names(nor_mat_BA3) = taiji_res$X
colnames(nor_mat_BA3) = c('GroupB3','GroupA3')
nor_mat_BA3$diff = nor_mat_BA3$GroupB3 - nor_mat_BA3$GroupA3
fit_diff_BA3 = fitdistr(abs(nor_mat_BA3$diff),'normal')
para_diff_BA3 = fit_diff_BA3$estimate
pvalue_diff_BA3 = pnorm(abs(nor_mat_BA3$diff), mean = para_diff_BA3[1], sd = para_diff_BA3[2],lower.tail = F)
FDR_diff_BA3 = p.adjust(pvalue_diff_BA3,method = "fdr", n = length(pvalue_diff_BA3))
nor_mat_BA3$pvalue = pvalue_diff_BA3
nor_mat_BA3$FDR = FDR_diff_BA3
nor_mat_BA3$TFs = row.names(nor_mat_BA3)
diff_TFs_BA3 = nor_mat_BA3[which(nor_mat_BA3$FDR < 0.05),]
new_matrix_BA3 = new_BA3[which(new_BA3$TFs %in% diff_TFs_BA3$TFs),]
final_diff_TFs_BA3 = inner_join(diff_TFs_BA3,new_matrix_BA3,by = 'TFs')
write.csv(final_diff_TFs_BA3,'Diff_TFs_BA3.csv',row.names = FALSE)

##################D1C#####################
new_D1C = bind_cols(taiji_res[,54],taiji_res[,56],taiji_res[,58],taiji_res[,60],taiji_res[,62],taiji_res[,64],
                    taiji_res[,42],taiji_res[,44],taiji_res[,46],taiji_res[,48],taiji_res[,50],taiji_res[,52])
colnames(new_D1C) = c('GroupD_1','GroupD_2','GroupD_3','GroupD_4','GroupD_5','GroupD_6',
                      'GroupC_1','GroupC_2','GroupC_3','GroupC_4','GroupC_5','GroupC_6')
new_D1C$TFs =taiji_res$X
row.names(new_D1C) = taiji_res$X
aver_GroupD1 = rowMeans(new_D1C[,1:6])
aver_GroupC1 = rowMeans(new_D1C[,7:12])
#####################potential key TFs#######################
fit_GroupD1 = fitdistr(aver_GroupD1,'normal')
para_GroupD1 = fit_GroupD1$estimate
pvalue_GroupD1 = pnorm(aver_GroupD1,mean = para_GroupD1[1], sd = para_GroupD1[2], lower.tail = F)
FDR_GroupD1 = p.adjust(pvalue_GroupD1,method = "fdr", n = length(pvalue_GroupD1))
p_GroupD1 = data.frame(aver_GroupD1,pvalue_GroupD1,FDR_GroupD1,taiji_res$X)
key_GroupD1 = p_GroupD1[which(p_GroupD1$FDR_GroupD1 < 0.05),]
write.csv(key_GroupD1,'potential_key_TFs_D1C_GroupD1.csv',row.names = FALSE)

fit_GroupC1 = fitdistr(aver_GroupC1,'normal')
para_GroupC1 = fit_GroupC1$estimate
pvalue_GroupC1 = pnorm(aver_GroupC1,mean = para_GroupC1[1], sd = para_GroupC1[2], lower.tail = F)
FDR_GroupC1 = p.adjust(pvalue_GroupC1,method = "fdr", n = length(pvalue_GroupC1))
p_GroupC1 = data.frame(aver_GroupC1,pvalue_GroupC1,FDR_GroupC1,taiji_res$X)
key_GroupC1 = p_GroupC1[which(p_GroupC1$FDR_GroupC1 < 0.05),]
write.csv(key_GroupC1,'potential_key_TFs_D1C_GroupC1.csv',row.names = FALSE)
#############################differential#######################
tmp = scale(data.frame(aver_GroupD1,aver_GroupC1))
nor_mat_D1C = data.frame(tmp)
row.names(nor_mat_D1C) = taiji_res$X
colnames(nor_mat_D1C) = c('GroupD1','GroupC1')
nor_mat_D1C$diff = nor_mat_D1C$GroupD1 - nor_mat_D1C$GroupC1
fit_diff_D1C = fitdistr(abs(nor_mat_D1C$diff),'normal')
para_diff_D1C = fit_diff_D1C$estimate
pvalue_diff_D1C = pnorm(abs(nor_mat_D1C$diff), mean = para_diff_D1C[1], sd = para_diff_D1C[2],lower.tail = F)
FDR_diff_D1C = p.adjust(pvalue_diff_D1C,method = "fdr", n = length(pvalue_diff_D1C))
nor_mat_D1C$pvalue = pvalue_diff_D1C
nor_mat_D1C$FDR = FDR_diff_D1C
nor_mat_D1C$TFs = row.names(nor_mat_D1C)
diff_TFs_D1C = nor_mat_D1C[which(nor_mat_D1C$FDR < 0.05),]
new_matrix_D1C = new_D1C[which(new_D1C$TFs %in% diff_TFs_D1C$TFs),]
final_diff_TFs_D1C = inner_join(diff_TFs_D1C,new_matrix_D1C,by = 'TFs')
write.csv(final_diff_TFs_D1C,'Diff_TFs_D1C.csv',row.names = FALSE)

##################D2C#####################
new_D2C = bind_cols(taiji_res[,66],taiji_res[,68],taiji_res[,70],
                    taiji_res[,42],taiji_res[,44],taiji_res[,46],taiji_res[,48],taiji_res[,50],taiji_res[,52])
colnames(new_D2C) = c('GroupD_7','GroupD_8','GroupD_9',
                      'GroupC_1','GroupC_2','GroupC_3','GroupC_4','GroupC_5','GroupC_6')
new_D2C$TFs =taiji_res$X
row.names(new_D2C) = taiji_res$X
aver_GroupD2 = rowMeans(new_D2C[,1:3])
aver_GroupC1 = rowMeans(new_D2C[,4:9])
#####################potential key TFs#######################
fit_GroupD2 = fitdistr(aver_GroupD2,'normal')
para_GroupD2 = fit_GroupD2$estimate
pvalue_GroupD2 = pnorm(aver_GroupD2,mean = para_GroupD2[1], sd = para_GroupD2[2], lower.tail = F)
FDR_GroupD2 = p.adjust(pvalue_GroupD2,method = "fdr", n = length(pvalue_GroupD2))
p_GroupD2 = data.frame(aver_GroupD2,pvalue_GroupD2,FDR_GroupD2,taiji_res$X)
key_GroupD2 = p_GroupD2[which(p_GroupD2$FDR_GroupD2 < 0.05),]
write.csv(key_GroupD2,'potential_key_TFs_D2C_GroupD2.csv',row.names = FALSE)

fit_GroupC1 = fitdistr(aver_GroupC1,'normal')
para_GroupC1 = fit_GroupC1$estimate
pvalue_GroupC1 = pnorm(aver_GroupC1,mean = para_GroupC1[1], sd = para_GroupC1[2], lower.tail = F)
FDR_GroupC1 = p.adjust(pvalue_GroupC1,method = "fdr", n = length(pvalue_GroupC1))
p_GroupC1 = data.frame(aver_GroupC1,pvalue_GroupC1,FDR_GroupC1,taiji_res$X)
key_GroupC1 = p_GroupC1[which(p_GroupC1$FDR_GroupC1 < 0.05),]
write.csv(key_GroupC1,'potential_key_TFs_D2C_GroupC1.csv',row.names = FALSE)
#############################differential#######################
tmp = scale(data.frame(aver_GroupD2,aver_GroupC1))
nor_mat_D2C = data.frame(tmp)
row.names(nor_mat_D2C) = taiji_res$X
colnames(nor_mat_D2C) = c('GroupD2','GroupC1')
nor_mat_D2C$diff = nor_mat_D2C$GroupD2 - nor_mat_D2C$GroupC1
fit_diff_D2C = fitdistr(abs(nor_mat_D2C$diff),'normal')
para_diff_D2C = fit_diff_D2C$estimate
pvalue_diff_D2C = pnorm(abs(nor_mat_D2C$diff), mean = para_diff_D2C[1], sd = para_diff_D2C[2],lower.tail = F)
FDR_diff_D2C = p.adjust(pvalue_diff_D2C,method = "fdr", n = length(pvalue_diff_D2C))
nor_mat_D2C$pvalue = pvalue_diff_D2C
nor_mat_D2C$FDR = FDR_diff_D2C
nor_mat_D2C$TFs = row.names(nor_mat_D2C)
diff_TFs_D2C = nor_mat_D2C[which(nor_mat_D2C$FDR < 0.05),]
new_matrix_D2C = new_D2C[which(new_D2C$TFs %in% diff_TFs_D2C$TFs),]
final_diff_TFs_D2C = inner_join(diff_TFs_D2C,new_matrix_D2C,by = 'TFs')
write.csv(final_diff_TFs_D2C,'Diff_TFs_D2C.csv',row.names = FALSE)

##################D3C#####################
new_D3C = bind_cols(taiji_res[,72],taiji_res[,74],taiji_res[,76],taiji_res[,77],taiji_res[,79],
                    taiji_res[,42],taiji_res[,44],taiji_res[,46],taiji_res[,48],taiji_res[,50],taiji_res[,52])
colnames(new_D3C) = c('GroupD_11','GroupD_12','GroupD_13','GroupD_14','GroupD_15',
                      'GroupC_1','GroupC_2','GroupC_3','GroupC_4','GroupC_5','GroupC_6')
new_D3C$TFs =taiji_res$X
row.names(new_D3C) = taiji_res$X
aver_GroupD3 = rowMeans(new_D3C[,1:5])
aver_GroupC1 = rowMeans(new_D3C[,6:11])
#####################potential key TFs#######################
fit_GroupD3 = fitdistr(aver_GroupD3,'normal')
para_GroupD3 = fit_GroupD3$estimate
pvalue_GroupD3 = pnorm(aver_GroupD3,mean = para_GroupD3[1], sd = para_GroupD3[2], lower.tail = F)
FDR_GroupD3 = p.adjust(pvalue_GroupD3,method = "fdr", n = length(pvalue_GroupD3))
p_GroupD3 = data.frame(aver_GroupD3,pvalue_GroupD3,FDR_GroupD3,taiji_res$X)
key_GroupD3 = p_GroupD3[which(p_GroupD3$FDR_GroupD3 < 0.05),]
write.csv(key_GroupD3,'potential_key_TFs_D3C_GroupD3.csv',row.names = FALSE)

fit_GroupC1 = fitdistr(aver_GroupC1,'normal')
para_GroupC1 = fit_GroupC1$estimate
pvalue_GroupC1 = pnorm(aver_GroupC1,mean = para_GroupC1[1], sd = para_GroupC1[2], lower.tail = F)
FDR_GroupC1 = p.adjust(pvalue_GroupC1,method = "fdr", n = length(pvalue_GroupC1))
p_GroupC1 = data.frame(aver_GroupC1,pvalue_GroupC1,FDR_GroupC1,taiji_res$X)
key_GroupC1 = p_GroupC1[which(p_GroupC1$FDR_GroupC1 < 0.05),]
write.csv(key_GroupC1,'potential_key_TFs_D3C_GroupC1.csv',row.names = FALSE)
#############################differential#######################
tmp = scale(data.frame(aver_GroupD3,aver_GroupC1))
nor_mat_D3C = data.frame(tmp)
row.names(nor_mat_D3C) = taiji_res$X
colnames(nor_mat_D3C) = c('GroupD3','GroupC1')
nor_mat_D3C$diff = nor_mat_D3C$GroupD3 - nor_mat_D3C$GroupC1
fit_diff_D3C = fitdistr(abs(nor_mat_D3C$diff),'normal')
para_diff_D3C = fit_diff_D3C$estimate
pvalue_diff_D3C = pnorm(abs(nor_mat_D3C$diff), mean = para_diff_D3C[1], sd = para_diff_D3C[2],lower.tail = F)
FDR_diff_D3C = p.adjust(pvalue_diff_D3C,method = "fdr", n = length(pvalue_diff_D3C))
nor_mat_D3C$pvalue = pvalue_diff_D3C
nor_mat_D3C$FDR = FDR_diff_D3C
nor_mat_D3C$TFs = row.names(nor_mat_D3C)
diff_TFs_D3C = nor_mat_D3C[which(nor_mat_D3C$FDR < 0.05),]
new_matrix_D3C = new_D3C[which(new_D3C$TFs %in% diff_TFs_D3C$TFs),]
final_diff_TFs_D3C = inner_join(diff_TFs_D3C,new_matrix_D3C,by = 'TFs')
write.csv(final_diff_TFs_D3C,'Diff_TFs_D3C.csv',row.names = FALSE)

##################CA1#####################
new_CA1 = bind_cols(taiji_res[,42],taiji_res[,44],taiji_res[,46],taiji_res[,48],taiji_res[,50],taiji_res[,52],
                     taiji_res[,2],taiji_res[,4],taiji_res[,6],taiji_res[,7],taiji_res[,9],taiji_res[,11])
colnames(new_CA1) = c('GroupC_1','GroupC_2','GroupC_3','GroupC_4','GroupC_5','GroupC_6',
                      'GroupA_1','GroupA_2','GroupA_3','GroupA_4','GroupA_5','GroupA_6')
new_CA1$TFs =taiji_res$X
row.names(new_CA1) = taiji_res$X
aver_GroupC = rowMeans(new_CA1[,1:6])
aver_GroupA1 = rowMeans(new_CA1[,7:12])
#####################potential key TFs#######################
fit_GroupC = fitdistr(aver_GroupC,'normal')
para_GroupC = fit_GroupC$estimate
pvalue_GroupC = pnorm(aver_GroupC,mean = para_GroupC[1], sd = para_GroupC[2], lower.tail = F)
FDR_GroupC = p.adjust(pvalue_GroupC,method = "fdr", n = length(pvalue_GroupC))
p_GroupC = data.frame(aver_GroupC,pvalue_GroupC,FDR_GroupC,taiji_res$X)
key_GroupC = p_GroupC[which(p_GroupC$FDR_GroupC < 0.05),]
write.csv(key_GroupC,'potential_key_TFs_CA1_GroupC.csv',row.names = FALSE)

fit_GroupA1 = fitdistr(aver_GroupA1,'normal')
para_GroupA1 = fit_GroupA1$estimate
pvalue_GroupA1 = pnorm(aver_GroupA1,mean = para_GroupA1[1], sd = para_GroupA1[2], lower.tail = F)
FDR_GroupA1 = p.adjust(pvalue_GroupA1,method = "fdr", n = length(pvalue_GroupA1))
p_GroupA1 = data.frame(aver_GroupA1,pvalue_GroupA1,FDR_GroupA1,taiji_res$X)
key_GroupA1 = p_GroupA1[which(p_GroupA1$FDR_GroupA1 < 0.05),]
write.csv(key_GroupA1,'potential_key_TFs_CA1_GroupA1.csv',row.names = FALSE)
#############################differential#######################
tmp = scale(data.frame(aver_GroupC,aver_GroupA1))
nor_mat_CA1 = data.frame(tmp)
row.names(nor_mat_CA1) = taiji_res$X
colnames(nor_mat_CA1) = c('GroupC','GroupA1')
nor_mat_CA1$diff = nor_mat_CA1$GroupC - nor_mat_CA1$GroupA1
fit_diff_CA1 = fitdistr(abs(nor_mat_CA1$diff),'normal')
para_diff_CA1 = fit_diff_CA1$estimate
pvalue_diff_CA1 = pnorm(abs(nor_mat_CA1$diff), mean = para_diff_CA1[1], sd = para_diff_CA1[2],lower.tail = F)
FDR_diff_CA1 = p.adjust(pvalue_diff_CA1,method = "fdr", n = length(pvalue_diff_CA1))
nor_mat_CA1$pvalue = pvalue_diff_CA1
nor_mat_CA1$FDR = FDR_diff_CA1
nor_mat_CA1$TFs = row.names(nor_mat_CA1)
diff_TFs_CA1 = nor_mat_CA1[which(nor_mat_CA1$FDR < 0.05),]
new_matrix_CA1 = new_CA1[which(new_CA1$TFs %in% diff_TFs_CA1$TFs),]
final_diff_TFs_CA1 = inner_join(diff_TFs_CA1,new_matrix_CA1,by = 'TFs')
write.csv(final_diff_TFs_CA1,'Diff_TFs_CA1.csv',row.names = FALSE)

##################CA2#####################
new_CA2 = bind_cols(taiji_res[,42],taiji_res[,44],taiji_res[,46],taiji_res[,48],taiji_res[,50],taiji_res[,52],
                     taiji_res[,13],taiji_res[,15],taiji_res[,17])
colnames(new_CA2) = c('GroupC_1','GroupC_2','GroupC_3','GroupC_4','GroupC_5','GroupC_6',
                       'GroupA_7','GroupA_8','GroupA_9')
new_CA2$TFs =taiji_res$X
row.names(new_CA2) = taiji_res$X
aver_GroupC = rowMeans(new_CA2[,1:6])
aver_GroupA2 = rowMeans(new_CA2[,7:9])
#####################potential key TFs#######################
fit_GroupC = fitdistr(aver_GroupC,'normal')
para_GroupC = fit_GroupC$estimate
pvalue_GroupC = pnorm(aver_GroupC,mean = para_GroupC[1], sd = para_GroupC[2], lower.tail = F)
FDR_GroupC = p.adjust(pvalue_GroupC,method = "fdr", n = length(pvalue_GroupC))
p_GroupC = data.frame(aver_GroupC,pvalue_GroupC,FDR_GroupC,taiji_res$X)
key_GroupC = p_GroupC[which(p_GroupC$FDR_GroupC < 0.05),]
write.csv(key_GroupC,'potential_key_TFs_CA2_GroupC.csv',row.names = FALSE)

fit_GroupA2 = fitdistr(aver_GroupA2,'normal')
para_GroupA2 = fit_GroupA2$estimate
pvalue_GroupA2 = pnorm(aver_GroupA2,mean = para_GroupA2[1], sd = para_GroupA2[2], lower.tail = F)
FDR_GroupA2 = p.adjust(pvalue_GroupA2,method = "fdr", n = length(pvalue_GroupA2))
p_GroupA2 = data.frame(aver_GroupA2,pvalue_GroupA2,FDR_GroupA2,taiji_res$X)
key_GroupA2 = p_GroupA2[which(p_GroupA2$FDR_GroupA2 < 0.05),]
write.csv(key_GroupA2,'potential_key_TFs_CA2_GroupA2.csv',row.names = FALSE)
#############################differential#######################
tmp = scale(data.frame(aver_GroupC,aver_GroupA2))
nor_mat_CA2 = data.frame(tmp)
row.names(nor_mat_CA2) = taiji_res$X
colnames(nor_mat_CA2) = c('GroupC','GroupA2')
nor_mat_CA2$diff = nor_mat_CA2$GroupC - nor_mat_CA2$GroupA2
fit_diff_CA2 = fitdistr(abs(nor_mat_CA2$diff),'normal')
para_diff_CA2 = fit_diff_CA2$estimate
pvalue_diff_CA2 = pnorm(abs(nor_mat_CA2$diff), mean = para_diff_CA2[1], sd = para_diff_CA2[2],lower.tail = F)
FDR_diff_CA2 = p.adjust(pvalue_diff_CA2,method = "fdr", n = length(pvalue_diff_CA2))
nor_mat_CA2$pvalue = pvalue_diff_CA2
nor_mat_CA2$FDR = FDR_diff_CA2
nor_mat_CA2$TFs = row.names(nor_mat_CA2)
diff_TFs_CA2 = nor_mat_CA2[which(nor_mat_CA2$FDR < 0.05),]
new_matrix_CA2 = new_CA2[which(new_CA2$TFs %in% diff_TFs_CA2$TFs),]
final_diff_TFs_CA2 = inner_join(diff_TFs_CA2,new_matrix_CA2,by = 'TFs')
write.csv(final_diff_TFs_CA2,'Diff_TFs_CA2.csv',row.names = FALSE)

##################CA3#####################
new_CA3 = bind_cols(taiji_res[,42],taiji_res[,44],taiji_res[,46],taiji_res[,48],taiji_res[,50],taiji_res[,52],
                     taiji_res[,19],taiji_res[,20],taiji_res[,22],taiji_res[,24],taiji_res[,26],taiji_res[,28])
colnames(new_CA3) = c('GroupC_1','GroupC_2','GroupC_3','GroupC_4','GroupC_5','GroupC_6',
                       'GroupA_10','GroupA_11','GroupA_12','GroupA_13','GroupA_14','GroupA_15')
new_CA3$TFs =taiji_res$X
row.names(new_CA3) = taiji_res$X
aver_GroupC = rowMeans(new_CA3[,1:6])
aver_GroupA3 = rowMeans(new_CA3[,7:12])
#####################potential key TFs#######################
fit_GroupC = fitdistr(aver_GroupC,'normal')
para_GroupC = fit_GroupC$estimate
pvalue_GroupC = pnorm(aver_GroupC,mean = para_GroupC[1], sd = para_GroupC[2], lower.tail = F)
FDR_GroupC = p.adjust(pvalue_GroupC,method = "fdr", n = length(pvalue_GroupC))
p_GroupC = data.frame(aver_GroupC,pvalue_GroupC,FDR_GroupC,taiji_res$X)
key_GroupC = p_GroupC[which(p_GroupC$FDR_GroupC < 0.05),]
write.csv(key_GroupC,'potential_key_TFs_CA3_GroupC.csv',row.names = FALSE)

fit_GroupA3 = fitdistr(aver_GroupA3,'normal')
para_GroupA3 = fit_GroupA3$estimate
pvalue_GroupA3 = pnorm(aver_GroupA3,mean = para_GroupA3[1], sd = para_GroupA3[2], lower.tail = F)
FDR_GroupA3 = p.adjust(pvalue_GroupA3,method = "fdr", n = length(pvalue_GroupA3))
p_GroupA3 = data.frame(aver_GroupA3,pvalue_GroupA3,FDR_GroupA3,taiji_res$X)
key_GroupA3 = p_GroupA3[which(p_GroupA3$FDR_GroupA3 < 0.05),]
write.csv(key_GroupA3,'potential_key_TFs_CA3_GroupA3.csv',row.names = FALSE)
#############################differential#######################
tmp = scale(data.frame(aver_GroupC,aver_GroupA3))
nor_mat_CA3 = data.frame(tmp)
row.names(nor_mat_CA3) = taiji_res$X
colnames(nor_mat_CA3) = c('GroupC','GroupA3')
nor_mat_CA3$diff = nor_mat_CA3$GroupC - nor_mat_CA3$GroupA3
fit_diff_CA3 = fitdistr(abs(nor_mat_CA3$diff),'normal')
para_diff_CA3 = fit_diff_CA3$estimate
pvalue_diff_CA3 = pnorm(abs(nor_mat_CA3$diff), mean = para_diff_CA3[1], sd = para_diff_CA3[2],lower.tail = F)
FDR_diff_CA3 = p.adjust(pvalue_diff_CA3,method = "fdr", n = length(pvalue_diff_CA3))
nor_mat_CA3$pvalue = pvalue_diff_CA3
nor_mat_CA3$FDR = FDR_diff_CA3
nor_mat_CA3$TFs = row.names(nor_mat_CA3)
diff_TFs_CA3 = nor_mat_CA3[which(nor_mat_CA3$FDR < 0.05),]
new_matrix_CA3 = new_CA3[which(new_CA3$TFs %in% diff_TFs_CA3$TFs),]
final_diff_TFs_CA3 = inner_join(diff_TFs_CA3,new_matrix_CA3,by = 'TFs')
write.csv(final_diff_TFs_CA3,'Diff_TFs_CA3.csv',row.names = FALSE)

##################BD1#####################
new_BD1 = bind_cols(taiji_res[,30],taiji_res[,32],taiji_res[,34],taiji_res[,36],taiji_res[,38],taiji_res[,40],
                    taiji_res[,54],taiji_res[,56],taiji_res[,58],taiji_res[,60],taiji_res[,62],taiji_res[,64])
colnames(new_BD1) = c('GroupB_1','GroupB_2','GroupB_3','GroupB_4','GroupB_5','GroupB_6',
                      'GroupD_1','GroupD_2','GroupD_3','GroupD_4','GroupD_5','GroupD_6')
new_BD1$TFs =taiji_res$X
row.names(new_BD1) = taiji_res$X
aver_GroupB1 = rowMeans(new_BD1[,1:6])
aver_GroupD1 = rowMeans(new_BD1[,7:12])
#####################potential key TFs#######################
fit_GroupB1 = fitdistr(aver_GroupB1,'normal')
para_GroupB1 = fit_GroupB1$estimate
pvalue_GroupB1 = pnorm(aver_GroupB1,mean = para_GroupB1[1], sd = para_GroupB1[2], lower.tail = F)
FDR_GroupB1 = p.adjust(pvalue_GroupB1,method = "fdr", n = length(pvalue_GroupB1))
p_GroupB1 = data.frame(aver_GroupB1,pvalue_GroupB1,FDR_GroupB1,taiji_res$X)
key_GroupB1 = p_GroupB1[which(p_GroupB1$FDR_GroupB1 < 0.05),]
write.csv(key_GroupB1,'potential_key_TFs_BD1_GroupB1.csv',row.names = FALSE)

fit_GroupD1 = fitdistr(aver_GroupD1,'normal')
para_GroupD1 = fit_GroupD1$estimate
pvalue_GroupD1 = pnorm(aver_GroupD1,mean = para_GroupD1[1], sd = para_GroupD1[2], lower.tail = F)
FDR_GroupD1 = p.adjust(pvalue_GroupD1,method = "fdr", n = length(pvalue_GroupD1))
p_GroupD1 = data.frame(aver_GroupD1,pvalue_GroupD1,FDR_GroupD1,taiji_res$X)
key_GroupD1 = p_GroupD1[which(p_GroupD1$FDR_GroupD1 < 0.05),]
write.csv(key_GroupD1,'potential_key_TFs_BD1_GroupD1.csv',row.names = FALSE)
#############################differential#######################
tmp = scale(data.frame(aver_GroupB1,aver_GroupD1))
nor_mat_BD1 = data.frame(tmp)
row.names(nor_mat_BD1) = taiji_res$X
colnames(nor_mat_BD1) = c('GroupB1','GroupD1')
nor_mat_BD1$diff = nor_mat_BD1$GroupB1 - nor_mat_BD1$GroupD1
fit_diff_BD1 = fitdistr(abs(nor_mat_BD1$diff),'normal')
para_diff_BD1 = fit_diff_BD1$estimate
pvalue_diff_BD1 = pnorm(abs(nor_mat_BD1$diff), mean = para_diff_BD1[1], sd = para_diff_BD1[2],lower.tail = F)
FDR_diff_BD1 = p.adjust(pvalue_diff_BD1,method = "fdr", n = length(pvalue_diff_BD1))
nor_mat_BD1$pvalue = pvalue_diff_BD1
nor_mat_BD1$FDR = FDR_diff_BD1
nor_mat_BD1$TFs = row.names(nor_mat_BD1)
diff_TFs_BD1 = nor_mat_BD1[which(nor_mat_BD1$FDR < 0.05),]
new_matrix_BD1 = new_BD1[which(new_BD1$TFs %in% diff_TFs_BD1$TFs),]
final_diff_TFs_BD1 = inner_join(diff_TFs_BD1,new_matrix_BD1,by = 'TFs')
write.csv(final_diff_TFs_BD1,'Diff_TFs_BD1.csv',row.names = FALSE)

##################BD2#####################
new_BD2 = bind_cols(taiji_res[,30],taiji_res[,32],taiji_res[,34],taiji_res[,36],taiji_res[,38],taiji_res[,40],
                    taiji_res[,66],taiji_res[,68],taiji_res[,70])
colnames(new_BD2) = c('GroupB_1','GroupB_2','GroupB_3','GroupB_4','GroupB_5','GroupB_6',
                      'GroupD_7','GroupD_8','GroupD_9')
new_BD2$TFs =taiji_res$X
row.names(new_BD2) = taiji_res$X
aver_GroupB2 = rowMeans(new_BD2[,1:6])
aver_GroupD2 = rowMeans(new_BD2[,7:9])
#####################potential key TFs#######################
fit_GroupB2 = fitdistr(aver_GroupB2,'normal')
para_GroupB2 = fit_GroupB2$estimate
pvalue_GroupB2 = pnorm(aver_GroupB2,mean = para_GroupB2[1], sd = para_GroupB2[2], lower.tail = F)
FDR_GroupB2 = p.adjust(pvalue_GroupB2,method = "fdr", n = length(pvalue_GroupB2))
p_GroupB2 = data.frame(aver_GroupB2,pvalue_GroupB2,FDR_GroupB2,taiji_res$X)
key_GroupB2 = p_GroupB2[which(p_GroupB2$FDR_GroupB2 < 0.05),]
write.csv(key_GroupB2,'potential_key_TFs_BD2_GroupB2.csv',row.names = FALSE)

fit_GroupD2 = fitdistr(aver_GroupD2,'normal')
para_GroupD2 = fit_GroupD2$estimate
pvalue_GroupD2 = pnorm(aver_GroupD2,mean = para_GroupD2[1], sd = para_GroupD2[2], lower.tail = F)
FDR_GroupD2 = p.adjust(pvalue_GroupD2,method = "fdr", n = length(pvalue_GroupD2))
p_GroupD2 = data.frame(aver_GroupD2,pvalue_GroupD2,FDR_GroupD2,taiji_res$X)
key_GroupD2 = p_GroupD2[which(p_GroupD2$FDR_GroupD2 < 0.05),]
write.csv(key_GroupD2,'potential_key_TFs_BD2_GroupD2.csv',row.names = FALSE)
#############################differential#######################
tmp = scale(data.frame(aver_GroupB2,aver_GroupD2))
nor_mat_BD2 = data.frame(tmp)
row.names(nor_mat_BD2) = taiji_res$X
colnames(nor_mat_BD2) = c('GroupB2','GroupD2')
nor_mat_BD2$diff = nor_mat_BD2$GroupB2 - nor_mat_BD2$GroupD2
fit_diff_BD2 = fitdistr(abs(nor_mat_BD2$diff),'normal')
para_diff_BD2 = fit_diff_BD2$estimate
pvalue_diff_BD2 = pnorm(abs(nor_mat_BD2$diff), mean = para_diff_BD2[1], sd = para_diff_BD2[2],lower.tail = F)
FDR_diff_BD2 = p.adjust(pvalue_diff_BD2,method = "fdr", n = length(pvalue_diff_BD2))
nor_mat_BD2$pvalue = pvalue_diff_BD2
nor_mat_BD2$FDR = FDR_diff_BD2
nor_mat_BD2$TFs = row.names(nor_mat_BD2)
diff_TFs_BD2 = nor_mat_BD2[which(nor_mat_BD2$FDR < 0.05),]
new_matrix_BD2 = new_BD2[which(new_BD2$TFs %in% diff_TFs_BD2$TFs),]
final_diff_TFs_BD2 = inner_join(diff_TFs_BD2,new_matrix_BD2,by = 'TFs')
write.csv(final_diff_TFs_BD2,'Diff_TFs_BD2.csv',row.names = FALSE)

##################BD3#####################
new_BD3 = bind_cols(taiji_res[,30],taiji_res[,32],taiji_res[,34],taiji_res[,36],taiji_res[,38],taiji_res[,40],
                    taiji_res[,72],taiji_res[,74],taiji_res[,76],taiji_res[,77],taiji_res[,79])
colnames(new_BD3) = c('GroupB_1','GroupB_2','GroupB_3','GroupB_4','GroupB_5','GroupB_6',
                      'GroupD_11','GroupD_12','GroupD_13','GroupD_14','GroupD_15')
new_BD3$TFs =taiji_res$X
row.names(new_BD3) = taiji_res$X
aver_GroupB3 = rowMeans(new_BD3[,1:6])
aver_GroupD3 = rowMeans(new_BD3[,7:11])
#####################potential key TFs#######################
fit_GroupB3 = fitdistr(aver_GroupB3,'normal')
para_GroupB3 = fit_GroupB3$estimate
pvalue_GroupB3 = pnorm(aver_GroupB3,mean = para_GroupB3[1], sd = para_GroupB3[2], lower.tail = F)
FDR_GroupB3 = p.adjust(pvalue_GroupB3,method = "fdr", n = length(pvalue_GroupB3))
p_GroupB3 = data.frame(aver_GroupB3,pvalue_GroupB3,FDR_GroupB3,taiji_res$X)
key_GroupB3 = p_GroupB3[which(p_GroupB3$FDR_GroupB3 < 0.05),]
write.csv(key_GroupB3,'potential_key_TFs_BD3_GroupB3.csv',row.names = FALSE)

fit_GroupD3 = fitdistr(aver_GroupD3,'normal')
para_GroupD3 = fit_GroupD3$estimate
pvalue_GroupD3 = pnorm(aver_GroupD3,mean = para_GroupD3[1], sd = para_GroupD3[2], lower.tail = F)
FDR_GroupD3 = p.adjust(pvalue_GroupD3,method = "fdr", n = length(pvalue_GroupD3))
p_GroupD3 = data.frame(aver_GroupD3,pvalue_GroupD3,FDR_GroupD3,taiji_res$X)
key_GroupD3 = p_GroupD3[which(p_GroupD3$FDR_GroupD3 < 0.05),]
write.csv(key_GroupD3,'potential_key_TFs_BD3_GroupD3.csv',row.names = FALSE)
#############################differential#######################
tmp = scale(data.frame(aver_GroupB3,aver_GroupD3))
nor_mat_BD3 = data.frame(tmp)
row.names(nor_mat_BD3) = taiji_res$X
colnames(nor_mat_BD3) = c('GroupB3','GroupD3')
nor_mat_BD3$diff = nor_mat_BD3$GroupB3 - nor_mat_BD3$GroupD3
fit_diff_BD3 = fitdistr(abs(nor_mat_BD3$diff),'normal')
para_diff_BD3 = fit_diff_BD3$estimate
pvalue_diff_BD3 = pnorm(abs(nor_mat_BD3$diff), mean = para_diff_BD3[1], sd = para_diff_BD3[2],lower.tail = F)
FDR_diff_BD3 = p.adjust(pvalue_diff_BD3,method = "fdr", n = length(pvalue_diff_BD3))
nor_mat_BD3$pvalue = pvalue_diff_BD3
nor_mat_BD3$FDR = FDR_diff_BD3
nor_mat_BD3$TFs = row.names(nor_mat_BD3)
diff_TFs_BD3 = nor_mat_BD3[which(nor_mat_BD3$FDR < 0.05),]
new_matrix_BD3 = new_BD3[which(new_BD3$TFs %in% diff_TFs_BD3$TFs),]
final_diff_TFs_BD3 = inner_join(diff_TFs_BD3,new_matrix_BD3,by = 'TFs')
write.csv(final_diff_TFs_BD3,'Diff_TFs_BD3.csv',row.names = FALSE)

##################A1E#####################

new_A1E = bind_cols(taiji_res[,2],taiji_res[,4],taiji_res[,6],taiji_res[,7],taiji_res[,9],taiji_res[,11],
                    taiji_res_E[,2],taiji_res_E[,4],taiji_res_E[,6],taiji_res_E[,8])
colnames(new_A1E) = c('GroupA_1','GroupA_2','GroupA_3','GroupA_4','GroupA_5','GroupA_6',
                      'GroupE_2','GroupE_4','GroupE_5','GroupE_6')
new_A1E$TFs =taiji_res$X
row.names(new_A1E) = taiji_res$X
aver_GroupA1 = rowMeans(new_A1E[,1:6])
aver_GroupE = rowMeans(new_A1E[,7:10])
#####################potential key TFs#######################
fit_GroupA1 = fitdistr(aver_GroupA1,'normal')
para_GroupA1 = fit_GroupA1$estimate
pvalue_GroupA1 = pnorm(aver_GroupA1,mean = para_GroupA1[1], sd = para_GroupA1[2], lower.tail = F)
FDR_GroupA1 = p.adjust(pvalue_GroupA1,method = "fdr", n = length(pvalue_GroupA1))
p_GroupA1 = data.frame(aver_GroupA1,pvalue_GroupA1,FDR_GroupA1,taiji_res$X)
key_GroupA1 = p_GroupA1[which(p_GroupA1$FDR_GroupA1 < 0.05),]
write.csv(key_GroupA1,'potential_key_TFs_A1E_GroupA1.csv',row.names = FALSE)

fit_GroupE = fitdistr(aver_GroupE,'normal')
para_GroupE = fit_GroupE$estimate
pvalue_GroupE = pnorm(aver_GroupE,mean = para_GroupE[1], sd = para_GroupE[2], lower.tail = F)
FDR_GroupE = p.adjust(pvalue_GroupE,method = "fdr", n = length(pvalue_GroupE))
p_GroupE = data.frame(aver_GroupE,pvalue_GroupE,FDR_GroupE,taiji_res$X)
key_GroupE = p_GroupE[which(p_GroupE$FDR_GroupE < 0.05),]
write.csv(key_GroupE,'potential_key_TFs_A1E_GroupE.csv',row.names = FALSE)
#############################differential#######################
tmp = scale(data.frame(aver_GroupA1,aver_GroupE))
nor_mat_A1E = data.frame(tmp)
row.names(nor_mat_A1E) = taiji_res$X
colnames(nor_mat_A1E) = c('GroupA1','GroupE')
nor_mat_A1E$diff = nor_mat_A1E$GroupA1 - nor_mat_A1E$GroupE
fit_diff_A1E = fitdistr(abs(nor_mat_A1E$diff),'normal')
para_diff_A1E = fit_diff_A1E$estimate
pvalue_diff_A1E = pnorm(abs(nor_mat_A1E$diff), mean = para_diff_A1E[1], sd = para_diff_A1E[2],lower.tail = F)
FDR_diff_A1E = p.adjust(pvalue_diff_A1E,method = "fdr", n = length(pvalue_diff_A1E))
nor_mat_A1E$pvalue = pvalue_diff_A1E
nor_mat_A1E$FDR = FDR_diff_A1E
nor_mat_A1E$TFs = row.names(nor_mat_A1E)
diff_TFs_A1E = nor_mat_A1E[which(nor_mat_A1E$FDR < 0.05),]
new_matrix_A1E = new_A1E[which(new_A1E$TFs %in% diff_TFs_A1E$TFs),]
final_diff_TFs_A1E = inner_join(diff_TFs_A1E,new_matrix_A1E,by = 'TFs')
write.csv(final_diff_TFs_A1E,'Diff_TFs_A1E.csv',row.names = FALSE)

##################A2E#####################
new_A2E = bind_cols(taiji_res[,13],taiji_res[,15],taiji_res[,17], taiji_res_E[,2],taiji_res_E[,4],taiji_res_E[,6],taiji_res_E[,8])
colnames(new_A2E) = c('GroupA_7','GroupA_8','GroupA_9','GroupE_2','GroupE_4','GroupE_5','GroupE_6')
new_A2E$TFs =taiji_res$X
row.names(new_A2E) = taiji_res$X
aver_GroupA2 = rowMeans(new_A2E[,1:3])
aver_GroupE = rowMeans(new_A2E[,4:7])
#####################potential key TFs#######################
fit_GroupA2 = fitdistr(aver_GroupA2,'normal')
para_GroupA2 = fit_GroupA2$estimate
pvalue_GroupA2 = pnorm(aver_GroupA2,mean = para_GroupA2[1], sd = para_GroupA2[2], lower.tail = F)
FDR_GroupA2 = p.adjust(pvalue_GroupA2,method = "fdr", n = length(pvalue_GroupA2))
p_GroupA2 = data.frame(aver_GroupA2,pvalue_GroupA2,FDR_GroupA2,taiji_res$X)
key_GroupA2 = p_GroupA2[which(p_GroupA2$FDR_GroupA2 < 0.05),]
write.csv(key_GroupA2,'potential_key_TFs_A2E_GroupA2.csv',row.names = FALSE)

fit_GroupE = fitdistr(aver_GroupE,'normal')
para_GroupE = fit_GroupE$estimate
pvalue_GroupE = pnorm(aver_GroupE,mean = para_GroupE[1], sd = para_GroupE[2], lower.tail = F)
FDR_GroupE = p.adjust(pvalue_GroupE,method = "fdr", n = length(pvalue_GroupE))
p_GroupE = data.frame(aver_GroupE,pvalue_GroupE,FDR_GroupE,taiji_res$X)
key_GroupE = p_GroupE[which(p_GroupE$FDR_GroupE < 0.05),]
write.csv(key_GroupE,'potential_key_TFs_A2E_GroupE.csv',row.names = FALSE)
#############################differential#######################
tmp = scale(data.frame(aver_GroupA2,aver_GroupE))
nor_mat_A2E = data.frame(tmp)
row.names(nor_mat_A2E) = taiji_res$X
colnames(nor_mat_A2E) = c('GroupA2','GroupE')
nor_mat_A2E$diff = nor_mat_A2E$GroupA2 - nor_mat_A2E$GroupE
fit_diff_A2E = fitdistr(abs(nor_mat_A2E$diff),'normal')
para_diff_A2E = fit_diff_A2E$estimate
pvalue_diff_A2E = pnorm(abs(nor_mat_A2E$diff), mean = para_diff_A2E[1], sd = para_diff_A2E[2],lower.tail = F)
FDR_diff_A2E = p.adjust(pvalue_diff_A2E,method = "fdr", n = length(pvalue_diff_A2E))
nor_mat_A2E$pvalue = pvalue_diff_A2E
nor_mat_A2E$FDR = FDR_diff_A2E
nor_mat_A2E$TFs = row.names(nor_mat_A2E)
diff_TFs_A2E = nor_mat_A2E[which(nor_mat_A2E$FDR < 0.05),]
new_matrix_A2E = new_A2E[which(new_A2E$TFs %in% diff_TFs_A2E$TFs),]
final_diff_TFs_A2E = inner_join(diff_TFs_A2E,new_matrix_A2E,by = 'TFs')
write.csv(final_diff_TFs_A2E,'Diff_TFs_A2E.csv',row.names = FALSE)

##################A3E#####################
new_A3E = bind_cols(taiji_res[,19],taiji_res[,20],taiji_res[,22],taiji_res[,24],
                    taiji_res[,26],taiji_res[,28], taiji_res_E[,2],taiji_res_E[,4],taiji_res_E[,6],taiji_res_E[,8])
colnames(new_A3E) = c('GroupA_10','GroupA_11','GroupA_12','GroupA_13',
                      'GroupA_14','GroupA_15','GroupE_2','GroupE_4','GroupE_5','GroupE_6')
new_A3E$TFs =taiji_res$X
row.names(new_A3E) = taiji_res$X
aver_GroupA3 = rowMeans(new_A3E[,1:6])
aver_GroupE = rowMeans(new_A3E[,7:10])
#####################potential key TFs#######################
fit_GroupA3 = fitdistr(aver_GroupA3,'normal')
para_GroupA3 = fit_GroupA3$estimate
pvalue_GroupA3 = pnorm(aver_GroupA3,mean = para_GroupA3[1], sd = para_GroupA3[2], lower.tail = F)
FDR_GroupA3 = p.adjust(pvalue_GroupA3,method = "fdr", n = length(pvalue_GroupA3))
p_GroupA3 = data.frame(aver_GroupA3,pvalue_GroupA3,FDR_GroupA3,taiji_res$X)
key_GroupA3 = p_GroupA3[which(p_GroupA3$FDR_GroupA3 < 0.05),]
write.csv(key_GroupA3,'potential_key_TFs_A3E_GroupA3.csv',row.names = FALSE)

fit_GroupE = fitdistr(aver_GroupE,'normal')
para_GroupE = fit_GroupE$estimate
pvalue_GroupE = pnorm(aver_GroupE,mean = para_GroupE[1], sd = para_GroupE[2], lower.tail = F)
FDR_GroupE = p.adjust(pvalue_GroupE,method = "fdr", n = length(pvalue_GroupE))
p_GroupE = data.frame(aver_GroupE,pvalue_GroupE,FDR_GroupE,taiji_res$X)
key_GroupE = p_GroupE[which(p_GroupE$FDR_GroupE < 0.05),]
write.csv(key_GroupE,'potential_key_TFs_A3E_GroupE.csv',row.names = FALSE)
#############################differential#######################
tmp = scale(data.frame(aver_GroupA3,aver_GroupE))
nor_mat_A3E = data.frame(tmp)
row.names(nor_mat_A3E) = taiji_res$X
colnames(nor_mat_A3E) = c('GroupA3','GroupE')
nor_mat_A3E$diff = nor_mat_A3E$GroupA3 - nor_mat_A3E$GroupE
fit_diff_A3E = fitdistr(abs(nor_mat_A3E$diff),'normal')
para_diff_A3E = fit_diff_A3E$estimate
pvalue_diff_A3E = pnorm(abs(nor_mat_A3E$diff), mean = para_diff_A3E[1], sd = para_diff_A3E[2],lower.tail = F)
FDR_diff_A3E = p.adjust(pvalue_diff_A3E,method = "fdr", n = length(pvalue_diff_A3E))
nor_mat_A3E$pvalue = pvalue_diff_A3E
nor_mat_A3E$FDR = FDR_diff_A3E
nor_mat_A3E$TFs = row.names(nor_mat_A3E)
diff_TFs_A3E = nor_mat_A3E[which(nor_mat_A3E$FDR < 0.05),]
new_matrix_A3E = new_A3E[which(new_A3E$TFs %in% diff_TFs_A3E$TFs),]
final_diff_TFs_A3E = inner_join(diff_TFs_A3E,new_matrix_A3E,by = 'TFs')
write.csv(final_diff_TFs_A3E,'Diff_TFs_A3E.csv',row.names = FALSE)

