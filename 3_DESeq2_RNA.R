###################set groups#############################
rm(list = ls())
options(stringsAsFactors = F)
library(dplyr)
library(DESeq2)
library(umap)
library(ggplot2)
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/RNA")
factors = read.csv('cofactors.csv')
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/RNA/count files")
Group_A = read.table('Group_A_id.txt',header = T)
rownames(Group_A) = Group_A[,1]
GroupA_1 = Group_A[,19:30]
GroupA_2 = Group_A[,31:36]
GroupA_3 = Group_A[,7:18]
Group_B = read.table('Group_B_id.txt',header = T)
rownames(Group_B) = Group_B[,1]
GroupB = Group_B[,7:18]
Group_C = read.table('Group_C_id.txt',header = T)
rownames(Group_C) = Group_C[,1]
GroupC = Group_C[,7:18]
Group_D = read.table('Group_D_id.txt',header = T)
rownames(Group_D) = Group_D[,1]
GroupD_1 = Group_D[,17:28]
GroupD_2 = Group_D[,29:34]
GroupD_3 = Group_D[,7:16]
Group_E = read.table('Group_E_id.txt',header = T)
rownames(Group_E) = Group_E[,1]
# GroupE = bind_cols(Group_E[,7:10],Group_E[,13:18])
GroupE = Group_E[,7:16]
Group_F = read.table('Group_F_id.txt',header = T)
rownames(Group_F) = Group_F[,1]
# GroupF = bind_cols(Group_F[,7:12],Group_F[,15:18])
GroupF = Group_F[,7:16]
Group_G = read.table('Group_G_id.txt',header = T)
rownames(Group_G) = Group_G[,1]
GroupG = Group_G[,7:18]

###############################B vs. A_1###############################################################
meta = bind_rows(factors[31:42,1:7],factors[1:12,1:7])
meta$condition_num = c(rep(1,12),rep(2,12))
group_BA1 = data.frame(GroupB,GroupA_1)
miss_1 <- c()
for (i in 1:nrow(group_BA1)) {
  if(length(which(group_BA1[i,1:12]<10)) > 0.5*ncol(group_BA1[,1:12])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(group_BA1)) {
  if(length(which(group_BA1[i,13:24]<10)) > 0.5*ncol(group_BA1[,13:24])) miss_2 <- append(miss_2,i)
}
miss = union(miss_1,miss_2)
group_BA1 = group_BA1[-miss,]
t_group_BA1 = t(group_BA1)
umap_results <- umap(t_group_BA1,n_neighbors=12,min_dist=0.1,n_components = 6)
umaps = data.frame(umap_results[["layout"]])
umaps$condition = as.factor(meta$condition)
ggplot(umaps,mapping = aes(x=X1,y=X2,col=condition)) +
  geom_point(aes(colour = condition))

cor_con_batch = cor(meta$condition_num,meta$batch,method = 'pearson')
cor_con_reads = cor(meta$condition_num,meta$reads,method = 'pearson')
cor_con_Align_rate = cor(meta$condition_num,meta$Align_rate,method = 'pearson')
cor_con_exon = cor(meta$condition_num,meta$exon,method = 'pearson')
cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
#  -con
# reads sequencing_depth -UMAP
corr_reads_1 = cor(umaps$X1,meta$reads,method = 'pearson')
corr_reads_2 = cor(umaps$X2,meta$reads,method = 'pearson')
corr_reads_3 = cor(umaps$X3,meta$reads,method = 'pearson')
corr_reads_4 = cor(umaps$X4,meta$reads,method = 'pearson')
corr_reads_5 = cor(umaps$X5,meta$reads,method = 'pearson')
corr_reads_6 = cor(umaps$X6,meta$reads,method = 'pearson')

corr_Align_rate_1 = cor(umaps$X1,meta$Align_rate,method = 'pearson')
corr_Align_rate_2 = cor(umaps$X2,meta$Align_rate,method = 'pearson')
corr_Align_rate_3 = cor(umaps$X3,meta$Align_rate,method = 'pearson')
corr_Align_rate_4 = cor(umaps$X4,meta$Align_rate,method = 'pearson')
corr_Align_rate_5 = cor(umaps$X5,meta$Align_rate,method = 'pearson')
corr_Align_rate_6 = cor(umaps$X6,meta$Align_rate,method = 'pearson')

corr_exon_1 = cor(umaps$X1,meta$exon,method = 'pearson')
corr_exon_2 = cor(umaps$X2,meta$exon,method = 'pearson')
corr_exon_3 = cor(umaps$X3,meta$exon,method = 'pearson')
corr_exon_4 = cor(umaps$X4,meta$exon,method = 'pearson')
corr_exon_5 = cor(umaps$X5,meta$exon,method = 'pearson')
corr_exon_6 = cor(umaps$X6,meta$exon,method = 'pearson')

corr_sequencing_depth_1 = cor(umaps$X1,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_2 = cor(umaps$X2,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_3 = cor(umaps$X3,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_4 = cor(umaps$X4,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_5 = cor(umaps$X5,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_6 = cor(umaps$X6,meta$sequencing_depth,method = 'pearson')

corr_batch_1 = cor(umaps$X1,meta$batch,method = 'pearson')
corr_batch_2 = cor(umaps$X2,meta$batch,method = 'pearson')
corr_batch_3 = cor(umaps$X3,meta$batch,method = 'pearson')
corr_batch_4 = cor(umaps$X4,meta$batch,method = 'pearson')
corr_batch_5 = cor(umaps$X5,meta$batch,method = 'pearson')
corr_batch_6 = cor(umaps$X6,meta$batch,method = 'pearson')

group_list1 = meta$condition
group_list2 = meta$reads
group_list3 = meta$Align_rate
group_list4 = meta$exon
group_list5 = meta$sequencing_depth
group_list6 = meta$batch

colData=data.frame(row.names = colnames(group_BA1),
                   condition=group_list1,
                   reads=group_list2,
                   Align_rate=group_list3,
                   exon = group_list4,
                   sequencing_depth = group_list5)

colData$condition = factor(colData$condition,levels = c('B','A'))
colData$reads =as.numeric(colData$reads)
colData$Align_rate=as.numeric(colData$Align_rate)
colData$exon=as.numeric(colData$exon)
colData$sequencing_depth=as.numeric(colData$sequencing_depth)

dds=DESeqDataSetFromMatrix(countData = group_BA1,
                           colData = colData,
                           design = ~ condition + reads + sequencing_depth)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","B","A"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.6,]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_BA1[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0),]
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/RNA/DEG")
all_glia = group_BA1[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_BA1.csv')
write.csv(DEG_matrix,'Differential_BA1.csv')
###############################B vs. A_2###############################################################
meta = bind_rows(factors[31:42,1:7],factors[13:18,1:7])
meta$condition_num = c(rep(1,12),rep(2,6))
group_BA2 = data.frame(GroupB,GroupA_2)
miss_1 <- c()
for (i in 1:nrow(group_BA2)) {
  if(length(which(group_BA2[i,1:12]<10)) > 0.5*ncol(group_BA2[,1:12])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(group_BA2)) {
  if(length(which(group_BA2[i,13:18]<10)) > 0.5*ncol(group_BA2[,13:18])) miss_2 <- append(miss_2,i)
}
miss = union(miss_1,miss_2)
group_BA2 = group_BA2[-miss,]
t_group_BA2 = t(group_BA2)
umap_results <- umap(t_group_BA2,n_neighbors=6,min_dist=0.1,n_components = 6)
umaps = data.frame(umap_results[["layout"]])
umaps$condition = as.factor(meta$condition)
ggplot(umaps,mapping = aes(x=X1,y=X2,col=condition)) +
  geom_point(aes(colour = condition))

cor_con_batch = cor(meta$condition_num,meta$batch,method = 'pearson')
cor_con_reads = cor(meta$condition_num,meta$reads,method = 'pearson')
cor_con_Align_rate = cor(meta$condition_num,meta$Align_rate,method = 'pearson')
cor_con_exon = cor(meta$condition_num,meta$exon,method = 'pearson')
cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
# exon -con
# exon sequencing_depth -UMAP
corr_reads_1 = cor(umaps$X1,meta$reads,method = 'pearson')
corr_reads_2 = cor(umaps$X2,meta$reads,method = 'pearson')
corr_reads_3 = cor(umaps$X3,meta$reads,method = 'pearson')
corr_reads_4 = cor(umaps$X4,meta$reads,method = 'pearson')
corr_reads_5 = cor(umaps$X5,meta$reads,method = 'pearson')
corr_reads_6 = cor(umaps$X6,meta$reads,method = 'pearson')

corr_Align_rate_1 = cor(umaps$X1,meta$Align_rate,method = 'pearson')
corr_Align_rate_2 = cor(umaps$X2,meta$Align_rate,method = 'pearson')
corr_Align_rate_3 = cor(umaps$X3,meta$Align_rate,method = 'pearson')
corr_Align_rate_4 = cor(umaps$X4,meta$Align_rate,method = 'pearson')
corr_Align_rate_5 = cor(umaps$X5,meta$Align_rate,method = 'pearson')
corr_Align_rate_6 = cor(umaps$X6,meta$Align_rate,method = 'pearson')

corr_exon_1 = cor(umaps$X1,meta$exon,method = 'pearson')
corr_exon_2 = cor(umaps$X2,meta$exon,method = 'pearson')
corr_exon_3 = cor(umaps$X3,meta$exon,method = 'pearson')
corr_exon_4 = cor(umaps$X4,meta$exon,method = 'pearson')
corr_exon_5 = cor(umaps$X5,meta$exon,method = 'pearson')
corr_exon_6 = cor(umaps$X6,meta$exon,method = 'pearson')

corr_sequencing_depth_1 = cor(umaps$X1,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_2 = cor(umaps$X2,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_3 = cor(umaps$X3,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_4 = cor(umaps$X4,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_5 = cor(umaps$X5,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_6 = cor(umaps$X6,meta$sequencing_depth,method = 'pearson')

corr_batch_1 = cor(umaps$X1,meta$batch,method = 'pearson')
corr_batch_2 = cor(umaps$X2,meta$batch,method = 'pearson')
corr_batch_3 = cor(umaps$X3,meta$batch,method = 'pearson')
corr_batch_4 = cor(umaps$X4,meta$batch,method = 'pearson')
corr_batch_5 = cor(umaps$X5,meta$batch,method = 'pearson')
corr_batch_6 = cor(umaps$X6,meta$batch,method = 'pearson')

group_list1 = meta$condition
group_list2 = meta$reads
group_list3 = meta$Align_rate
group_list4 = meta$exon
group_list5 = meta$sequencing_depth
group_list6 = meta$batch

colData=data.frame(row.names = colnames(group_BA2),
                   condition=group_list1,
                   reads=group_list2,
                   Align_rate=group_list3,
                   exon = group_list4,
                   sequencing_depth = group_list5)

colData$condition = factor(colData$condition,levels = c('B','A'))
colData$reads =as.numeric(colData$reads)
colData$Align_rate=as.numeric(colData$Align_rate)
colData$exon=as.numeric(colData$exon)
colData$sequencing_depth=as.numeric(colData$sequencing_depth)

dds=DESeqDataSetFromMatrix(countData = group_BA2,
                           colData = colData,
                           design = ~ condition + exon)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","B","A"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.6,]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_BA2[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0),]
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/RNA/DEG")
all_glia = group_BA2[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_BA2.csv')
write.csv(DEG_matrix,'Differential_BA2.csv')
###############################B vs. A_3###############################################################
meta = bind_rows(factors[31:42,1:7],factors[19:30,1:7])
meta$condition_num = c(rep(1,12),rep(2,12))
group_BA3 = data.frame(GroupB,GroupA_3)
miss_1 <- c()
for (i in 1:nrow(group_BA3)) {
  if(length(which(group_BA3[i,1:12]<10)) > 0.5*ncol(group_BA3[,1:12])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(group_BA3)) {
  if(length(which(group_BA3[i,13:24]<10)) > 0.5*ncol(group_BA3[,13:24])) miss_2 <- append(miss_2,i)
}
miss = union(miss_1,miss_2)
group_BA3 = group_BA3[-miss,]
t_group_BA3 = t(group_BA3)
umap_results <- umap(t_group_BA3,n_neighbors=12,min_dist=0.1,n_components = 6)
umaps = data.frame(umap_results[["layout"]])
umaps$condition = as.factor(meta$condition)
ggplot(umaps,mapping = aes(x=X1,y=X2,col=condition)) +
  geom_point(aes(colour = condition))

cor_con_batch = cor(meta$condition_num,meta$batch,method = 'pearson')
cor_con_reads = cor(meta$condition_num,meta$reads,method = 'pearson')
cor_con_Align_rate = cor(meta$condition_num,meta$Align_rate,method = 'pearson')
cor_con_exon = cor(meta$condition_num,meta$exon,method = 'pearson')
cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
# exon reads -con
# exon reads -UMAP
corr_reads_1 = cor(umaps$X1,meta$reads,method = 'pearson')
corr_reads_2 = cor(umaps$X2,meta$reads,method = 'pearson')
corr_reads_3 = cor(umaps$X3,meta$reads,method = 'pearson')
corr_reads_4 = cor(umaps$X4,meta$reads,method = 'pearson')
corr_reads_5 = cor(umaps$X5,meta$reads,method = 'pearson')
corr_reads_6 = cor(umaps$X6,meta$reads,method = 'pearson')

corr_Align_rate_1 = cor(umaps$X1,meta$Align_rate,method = 'pearson')
corr_Align_rate_2 = cor(umaps$X2,meta$Align_rate,method = 'pearson')
corr_Align_rate_3 = cor(umaps$X3,meta$Align_rate,method = 'pearson')
corr_Align_rate_4 = cor(umaps$X4,meta$Align_rate,method = 'pearson')
corr_Align_rate_5 = cor(umaps$X5,meta$Align_rate,method = 'pearson')
corr_Align_rate_6 = cor(umaps$X6,meta$Align_rate,method = 'pearson')

corr_exon_1 = cor(umaps$X1,meta$exon,method = 'pearson')
corr_exon_2 = cor(umaps$X2,meta$exon,method = 'pearson')
corr_exon_3 = cor(umaps$X3,meta$exon,method = 'pearson')
corr_exon_4 = cor(umaps$X4,meta$exon,method = 'pearson')
corr_exon_5 = cor(umaps$X5,meta$exon,method = 'pearson')
corr_exon_6 = cor(umaps$X6,meta$exon,method = 'pearson')

corr_sequencing_depth_1 = cor(umaps$X1,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_2 = cor(umaps$X2,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_3 = cor(umaps$X3,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_4 = cor(umaps$X4,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_5 = cor(umaps$X5,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_6 = cor(umaps$X6,meta$sequencing_depth,method = 'pearson')

corr_batch_1 = cor(umaps$X1,meta$batch,method = 'pearson')
corr_batch_2 = cor(umaps$X2,meta$batch,method = 'pearson')
corr_batch_3 = cor(umaps$X3,meta$batch,method = 'pearson')
corr_batch_4 = cor(umaps$X4,meta$batch,method = 'pearson')
corr_batch_5 = cor(umaps$X5,meta$batch,method = 'pearson')
corr_batch_6 = cor(umaps$X6,meta$batch,method = 'pearson')

group_list1 = meta$condition
group_list2 = meta$reads
group_list3 = meta$Align_rate
group_list4 = meta$exon
group_list5 = meta$sequencing_depth
group_list6 = meta$batch

colData=data.frame(row.names = colnames(group_BA3),
                   condition=group_list1,
                   reads=group_list2,
                   Align_rate=group_list3,
                   exon = group_list4,
                   sequencing_depth = group_list5)

colData$condition = factor(colData$condition,levels = c('B','A'))
colData$reads =as.numeric(colData$reads)
colData$Align_rate=as.numeric(colData$Align_rate)
colData$exon=as.numeric(colData$exon)
colData$sequencing_depth=as.numeric(colData$sequencing_depth)

dds=DESeqDataSetFromMatrix(countData = group_BA3,
                           colData = colData,
                           design = ~ condition + exon + reads)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","B","A"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.6,]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_BA3[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0),]
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/RNA/DEG")
all_glia = group_BA3[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_BA3.csv')
write.csv(DEG_matrix,'Differential_BA3.csv')
###############################D_1 vs. C###############################################################
meta = bind_rows(factors[55:66,1:7],factors[43:54,1:7])
meta$condition_num = c(rep(1,12),rep(2,12))
group_D1C = data.frame(GroupD_1,GroupC)
miss_1 <- c()
for (i in 1:nrow(group_D1C)) {
  if(length(which(group_D1C[i,1:12]<10)) > 0.5*ncol(group_D1C[,1:12])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(group_D1C)) {
  if(length(which(group_D1C[i,13:24]<10)) > 0.5*ncol(group_D1C[,13:24])) miss_2 <- append(miss_2,i)
}
miss = union(miss_1,miss_2)
group_D1C = group_D1C[-miss,]
t_group_D1C = t(group_D1C)
umap_results <- umap(t_group_D1C,n_neighbors=12,min_dist=0.1,n_components = 6)
umaps = data.frame(umap_results[["layout"]])
umaps$condition = as.factor(meta$condition)
ggplot(umaps,mapping = aes(x=X1,y=X2,col=condition)) +
  geom_point(aes(colour = condition))

cor_con_batch = cor(meta$condition_num,meta$batch,method = 'pearson')
cor_con_reads = cor(meta$condition_num,meta$reads,method = 'pearson')
cor_con_Align_rate = cor(meta$condition_num,meta$Align_rate,method = 'pearson')
cor_con_exon = cor(meta$condition_num,meta$exon,method = 'pearson')
cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
# sequencing_depth -con
# reads sequencing_depth -UMAP
corr_reads_1 = cor(umaps$X1,meta$reads,method = 'pearson')
corr_reads_2 = cor(umaps$X2,meta$reads,method = 'pearson')
corr_reads_3 = cor(umaps$X3,meta$reads,method = 'pearson')
corr_reads_4 = cor(umaps$X4,meta$reads,method = 'pearson')
corr_reads_5 = cor(umaps$X5,meta$reads,method = 'pearson')
corr_reads_6 = cor(umaps$X6,meta$reads,method = 'pearson')

corr_Align_rate_1 = cor(umaps$X1,meta$Align_rate,method = 'pearson')
corr_Align_rate_2 = cor(umaps$X2,meta$Align_rate,method = 'pearson')
corr_Align_rate_3 = cor(umaps$X3,meta$Align_rate,method = 'pearson')
corr_Align_rate_4 = cor(umaps$X4,meta$Align_rate,method = 'pearson')
corr_Align_rate_5 = cor(umaps$X5,meta$Align_rate,method = 'pearson')
corr_Align_rate_6 = cor(umaps$X6,meta$Align_rate,method = 'pearson')

corr_exon_1 = cor(umaps$X1,meta$exon,method = 'pearson')
corr_exon_2 = cor(umaps$X2,meta$exon,method = 'pearson')
corr_exon_3 = cor(umaps$X3,meta$exon,method = 'pearson')
corr_exon_4 = cor(umaps$X4,meta$exon,method = 'pearson')
corr_exon_5 = cor(umaps$X5,meta$exon,method = 'pearson')
corr_exon_6 = cor(umaps$X6,meta$exon,method = 'pearson')

corr_sequencing_depth_1 = cor(umaps$X1,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_2 = cor(umaps$X2,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_3 = cor(umaps$X3,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_4 = cor(umaps$X4,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_5 = cor(umaps$X5,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_6 = cor(umaps$X6,meta$sequencing_depth,method = 'pearson')

corr_batch_1 = cor(umaps$X1,meta$batch,method = 'pearson')
corr_batch_2 = cor(umaps$X2,meta$batch,method = 'pearson')
corr_batch_3 = cor(umaps$X3,meta$batch,method = 'pearson')
corr_batch_4 = cor(umaps$X4,meta$batch,method = 'pearson')
corr_batch_5 = cor(umaps$X5,meta$batch,method = 'pearson')
corr_batch_6 = cor(umaps$X6,meta$batch,method = 'pearson')

group_list1 = meta$condition
group_list2 = meta$reads
group_list3 = meta$Align_rate
group_list4 = meta$exon
group_list5 = meta$sequencing_depth
group_list6 = meta$batch

colData=data.frame(row.names = colnames(group_D1C),
                   condition=group_list1,
                   reads=group_list2,
                   Align_rate=group_list3,
                   exon = group_list4,
                   sequencing_depth = group_list5)

colData$condition = factor(colData$condition,levels = c('D','C'))
colData$reads =as.numeric(colData$reads)
colData$Align_rate=as.numeric(colData$Align_rate)
colData$exon=as.numeric(colData$exon)
colData$sequencing_depth=as.numeric(colData$sequencing_depth)

dds=DESeqDataSetFromMatrix(countData = group_D1C,
                           colData = colData,
                           design = ~ condition + reads + sequencing_depth)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","D","C"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.6,]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_D1C[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0),]

all_glia = group_D1C[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_D1C.csv')
write.csv(DEG_matrix,'Differential_D1C.csv')
###############################D_2 vs. C###############################################################
meta = bind_rows(factors[67:72,1:7],factors[43:54,1:7])
meta$condition_num = c(rep(1,6),rep(2,12))
group_D2C = data.frame(GroupD_2,GroupC)
miss_1 <- c()
for (i in 1:nrow(group_D2C)) {
  if(length(which(group_D2C[i,1:6]<10)) > 0.5*ncol(group_D2C[,1:6])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(group_D2C)) {
  if(length(which(group_D2C[i,7:18]<10)) > 0.5*ncol(group_D2C[,7:18])) miss_2 <- append(miss_2,i)
}
miss = union(miss_1,miss_2)
group_D2C = group_D2C[-miss,]
t_group_D2C = t(group_D2C)
umap_results <- umap(t_group_D2C,n_neighbors=6,min_dist=0.1,n_components = 6)
umaps = data.frame(umap_results[["layout"]])
umaps$condition = as.factor(meta$condition)
ggplot(umaps,mapping = aes(x=X1,y=X2,col=condition)) +
  geom_point(aes(colour = condition))

cor_con_batch = cor(meta$condition_num,meta$batch,method = 'pearson')
cor_con_reads = cor(meta$condition_num,meta$reads,method = 'pearson')
cor_con_Align_rate = cor(meta$condition_num,meta$Align_rate,method = 'pearson')
cor_con_exon = cor(meta$condition_num,meta$exon,method = 'pearson')
cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
# -con
# reads sequencing_depth -UMAP
corr_reads_1 = cor(umaps$X1,meta$reads,method = 'pearson')
corr_reads_2 = cor(umaps$X2,meta$reads,method = 'pearson')
corr_reads_3 = cor(umaps$X3,meta$reads,method = 'pearson')
corr_reads_4 = cor(umaps$X4,meta$reads,method = 'pearson')
corr_reads_5 = cor(umaps$X5,meta$reads,method = 'pearson')
corr_reads_6 = cor(umaps$X6,meta$reads,method = 'pearson')

corr_Align_rate_1 = cor(umaps$X1,meta$Align_rate,method = 'pearson')
corr_Align_rate_2 = cor(umaps$X2,meta$Align_rate,method = 'pearson')
corr_Align_rate_3 = cor(umaps$X3,meta$Align_rate,method = 'pearson')
corr_Align_rate_4 = cor(umaps$X4,meta$Align_rate,method = 'pearson')
corr_Align_rate_5 = cor(umaps$X5,meta$Align_rate,method = 'pearson')
corr_Align_rate_6 = cor(umaps$X6,meta$Align_rate,method = 'pearson')

corr_exon_1 = cor(umaps$X1,meta$exon,method = 'pearson')
corr_exon_2 = cor(umaps$X2,meta$exon,method = 'pearson')
corr_exon_3 = cor(umaps$X3,meta$exon,method = 'pearson')
corr_exon_4 = cor(umaps$X4,meta$exon,method = 'pearson')
corr_exon_5 = cor(umaps$X5,meta$exon,method = 'pearson')
corr_exon_6 = cor(umaps$X6,meta$exon,method = 'pearson')

corr_sequencing_depth_1 = cor(umaps$X1,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_2 = cor(umaps$X2,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_3 = cor(umaps$X3,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_4 = cor(umaps$X4,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_5 = cor(umaps$X5,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_6 = cor(umaps$X6,meta$sequencing_depth,method = 'pearson')

corr_batch_1 = cor(umaps$X1,meta$batch,method = 'pearson')
corr_batch_2 = cor(umaps$X2,meta$batch,method = 'pearson')
corr_batch_3 = cor(umaps$X3,meta$batch,method = 'pearson')
corr_batch_4 = cor(umaps$X4,meta$batch,method = 'pearson')
corr_batch_5 = cor(umaps$X5,meta$batch,method = 'pearson')
corr_batch_6 = cor(umaps$X6,meta$batch,method = 'pearson')

group_list1 = meta$condition
group_list2 = meta$reads
group_list3 = meta$Align_rate
group_list4 = meta$exon
group_list5 = meta$sequencing_depth
group_list6 = meta$batch

colData=data.frame(row.names = colnames(group_D2C),
                   condition=group_list1,
                   reads=group_list2,
                   Align_rate=group_list3,
                   exon = group_list4,
                   sequencing_depth = group_list5,
                   batch = group_list6)

colData$condition = factor(colData$condition,levels = c('D','C'))
colData$reads =as.numeric(colData$reads)
colData$Align_rate=as.numeric(colData$Align_rate)
colData$exon=as.numeric(colData$exon)
colData$sequencing_depth=as.numeric(colData$sequencing_depth)
colData$batch=as.numeric(colData$batch)

dds=DESeqDataSetFromMatrix(countData = group_D2C,
                           colData = colData,
                           design = ~ condition + reads)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","D","C"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.6,]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_D2C[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0),]

all_glia = group_D2C[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_D2C.csv')
write.csv(DEG_matrix,'Differential_D2C.csv')
###############################D_3 vs. C###############################################################
meta = bind_rows(factors[73:82,1:7],factors[43:54,1:7])
meta$condition_num = c(rep(1,10),rep(2,12))
group_D3C = data.frame(GroupD_3,GroupC)
miss_1 <- c()
for (i in 1:nrow(group_D3C)) {
  if(length(which(group_D3C[i,1:10]<10)) > 0.5*ncol(group_D3C[,1:10])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(group_D3C)) {
  if(length(which(group_D3C[i,11:22]<10)) > 0.5*ncol(group_D3C[,11:22])) miss_2 <- append(miss_2,i)
}
miss = union(miss_1,miss_2)
group_D3C = group_D3C[-miss,]
t_group_D3C = t(group_D3C)
umap_results <- umap(t_group_D3C,n_neighbors=12,min_dist=0.1,n_components = 6)
umaps = data.frame(umap_results[["layout"]])
umaps$condition = as.factor(meta$condition)
ggplot(umaps,mapping = aes(x=X1,y=X2,col=condition)) +
  geom_point(aes(colour = condition))

cor_con_batch = cor(meta$condition_num,meta$batch,method = 'pearson')
cor_con_reads = cor(meta$condition_num,meta$reads,method = 'pearson')
cor_con_Align_rate = cor(meta$condition_num,meta$Align_rate,method = 'pearson')
cor_con_exon = cor(meta$condition_num,meta$exon,method = 'pearson')
cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
# exon reads -con
# exon reads sequencing_depth -UMAP
corr_reads_1 = cor(umaps$X1,meta$reads,method = 'pearson')
corr_reads_2 = cor(umaps$X2,meta$reads,method = 'pearson')
corr_reads_3 = cor(umaps$X3,meta$reads,method = 'pearson')
corr_reads_4 = cor(umaps$X4,meta$reads,method = 'pearson')
corr_reads_5 = cor(umaps$X5,meta$reads,method = 'pearson')
corr_reads_6 = cor(umaps$X6,meta$reads,method = 'pearson')

corr_Align_rate_1 = cor(umaps$X1,meta$Align_rate,method = 'pearson')
corr_Align_rate_2 = cor(umaps$X2,meta$Align_rate,method = 'pearson')
corr_Align_rate_3 = cor(umaps$X3,meta$Align_rate,method = 'pearson')
corr_Align_rate_4 = cor(umaps$X4,meta$Align_rate,method = 'pearson')
corr_Align_rate_5 = cor(umaps$X5,meta$Align_rate,method = 'pearson')
corr_Align_rate_6 = cor(umaps$X6,meta$Align_rate,method = 'pearson')

corr_exon_1 = cor(umaps$X1,meta$exon,method = 'pearson')
corr_exon_2 = cor(umaps$X2,meta$exon,method = 'pearson')
corr_exon_3 = cor(umaps$X3,meta$exon,method = 'pearson')
corr_exon_4 = cor(umaps$X4,meta$exon,method = 'pearson')
corr_exon_5 = cor(umaps$X5,meta$exon,method = 'pearson')
corr_exon_6 = cor(umaps$X6,meta$exon,method = 'pearson')

corr_sequencing_depth_1 = cor(umaps$X1,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_2 = cor(umaps$X2,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_3 = cor(umaps$X3,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_4 = cor(umaps$X4,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_5 = cor(umaps$X5,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_6 = cor(umaps$X6,meta$sequencing_depth,method = 'pearson')

corr_batch_1 = cor(umaps$X1,meta$batch,method = 'pearson')
corr_batch_2 = cor(umaps$X2,meta$batch,method = 'pearson')
corr_batch_3 = cor(umaps$X3,meta$batch,method = 'pearson')
corr_batch_4 = cor(umaps$X4,meta$batch,method = 'pearson')
corr_batch_5 = cor(umaps$X5,meta$batch,method = 'pearson')
corr_batch_6 = cor(umaps$X6,meta$batch,method = 'pearson')

group_list1 = meta$condition
group_list2 = meta$reads
group_list3 = meta$Align_rate
group_list4 = meta$exon
group_list5 = meta$sequencing_depth
group_list6 = meta$batch

colData=data.frame(row.names = colnames(group_D3C),
                   condition=group_list1,
                   reads=group_list2,
                   Align_rate=group_list3,
                   exon = group_list4,
                   sequencing_depth = group_list5)

colData$condition = factor(colData$condition,levels = c('D','C'))
colData$reads =as.numeric(colData$reads)
colData$Align_rate=as.numeric(colData$Align_rate)
colData$exon=as.numeric(colData$exon)
colData$sequencing_depth=as.numeric(colData$sequencing_depth)

dds=DESeqDataSetFromMatrix(countData = group_D3C,
                           colData = colData,
                           design = ~ condition + reads + exon)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","D","C"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.6,]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_D3C[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0),]

all_glia = group_D3C[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_D3C.csv')
write.csv(DEG_matrix,'Differential_D3C.csv')
###############################C vs. A_1###############################################################
meta = bind_rows(factors[43:54,1:7],factors[1:12,1:7])
meta$condition_num = c(rep(1,12),rep(2,12))
group_CA1 = data.frame(GroupC,GroupA_1)
miss_1 <- c()
for (i in 1:nrow(group_CA1)) {
  if(length(which(group_CA1[i,1:12]<10)) > 0.5*ncol(group_CA1[,1:12])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(group_CA1)) {
  if(length(which(group_CA1[i,13:24]<10)) > 0.5*ncol(group_CA1[,13:24])) miss_2 <- append(miss_2,i)
}
miss = union(miss_1,miss_2)
group_CA1 = group_CA1[-miss,]
t_group_CA1 = t(group_CA1)
umap_results <- umap(t_group_CA1,n_neighbors=12,min_dist=0.1,n_components = 6)
umaps = data.frame(umap_results[["layout"]])
umaps$condition = as.factor(meta$condition)
ggplot(umaps,mapping = aes(x=X1,y=X2,col=condition)) +
  geom_point(aes(colour = condition))

cor_con_batch = cor(meta$condition_num,meta$batch,method = 'pearson')
cor_con_reads = cor(meta$condition_num,meta$reads,method = 'pearson')
cor_con_Align_rate = cor(meta$condition_num,meta$Align_rate,method = 'pearson')
cor_con_exon = cor(meta$condition_num,meta$exon,method = 'pearson')
cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
# -con
# reads sequencing_depth -UMAP
corr_reads_1 = cor(umaps$X1,meta$reads,method = 'pearson')
corr_reads_2 = cor(umaps$X2,meta$reads,method = 'pearson')
corr_reads_3 = cor(umaps$X3,meta$reads,method = 'pearson')
corr_reads_4 = cor(umaps$X4,meta$reads,method = 'pearson')
corr_reads_5 = cor(umaps$X5,meta$reads,method = 'pearson')
corr_reads_6 = cor(umaps$X6,meta$reads,method = 'pearson')

corr_Align_rate_1 = cor(umaps$X1,meta$Align_rate,method = 'pearson')
corr_Align_rate_2 = cor(umaps$X2,meta$Align_rate,method = 'pearson')
corr_Align_rate_3 = cor(umaps$X3,meta$Align_rate,method = 'pearson')
corr_Align_rate_4 = cor(umaps$X4,meta$Align_rate,method = 'pearson')
corr_Align_rate_5 = cor(umaps$X5,meta$Align_rate,method = 'pearson')
corr_Align_rate_6 = cor(umaps$X6,meta$Align_rate,method = 'pearson')

corr_exon_1 = cor(umaps$X1,meta$exon,method = 'pearson')
corr_exon_2 = cor(umaps$X2,meta$exon,method = 'pearson')
corr_exon_3 = cor(umaps$X3,meta$exon,method = 'pearson')
corr_exon_4 = cor(umaps$X4,meta$exon,method = 'pearson')
corr_exon_5 = cor(umaps$X5,meta$exon,method = 'pearson')
corr_exon_6 = cor(umaps$X6,meta$exon,method = 'pearson')

corr_sequencing_depth_1 = cor(umaps$X1,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_2 = cor(umaps$X2,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_3 = cor(umaps$X3,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_4 = cor(umaps$X4,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_5 = cor(umaps$X5,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_6 = cor(umaps$X6,meta$sequencing_depth,method = 'pearson')

corr_batch_1 = cor(umaps$X1,meta$batch,method = 'pearson')
corr_batch_2 = cor(umaps$X2,meta$batch,method = 'pearson')
corr_batch_3 = cor(umaps$X3,meta$batch,method = 'pearson')
corr_batch_4 = cor(umaps$X4,meta$batch,method = 'pearson')
corr_batch_5 = cor(umaps$X5,meta$batch,method = 'pearson')
corr_batch_6 = cor(umaps$X6,meta$batch,method = 'pearson')

group_list1 = meta$condition
group_list2 = meta$reads
group_list3 = meta$Align_rate
group_list4 = meta$exon
group_list5 = meta$sequencing_depth
group_list6 = meta$batch

colData=data.frame(row.names = colnames(group_CA1),
                   condition=group_list1,
                   reads=group_list2,
                   Align_rate=group_list3,
                   exon = group_list4,
                   sequencing_depth = group_list5,
                   batch = group_list6)

colData$condition = factor(colData$condition,levels = c('C','A'))
colData$reads =as.numeric(colData$reads)
colData$Align_rate=as.numeric(colData$Align_rate)
colData$exon=as.numeric(colData$exon)
colData$sequencing_depth=as.numeric(colData$sequencing_depth)
colData$batch=as.factor(colData$batch)

dds=DESeqDataSetFromMatrix(countData = group_CA1,
                           colData = colData,
                           design = ~ condition + reads)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","C","A"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.6,]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_CA1[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0),]
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/RNA/DEG")
all_glia = group_CA1[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_CA1.csv')
write.csv(DEG_matrix,'Differential_CA1.csv')
###############################C vs. A_2###############################################################
meta = bind_rows(factors[43:54,1:7],factors[13:18,1:7])
meta$condition_num = c(rep(1,12),rep(2,6))
group_CA2 = data.frame(GroupC,GroupA_2)
miss_1 <- c()
for (i in 1:nrow(group_CA2)) {
  if(length(which(group_CA2[i,1:12]<10)) > 0.5*ncol(group_CA2[,1:12])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(group_CA2)) {
  if(length(which(group_CA2[i,13:18]<10)) > 0.5*ncol(group_CA2[,13:18])) miss_2 <- append(miss_2,i)
}
miss = union(miss_1,miss_2)
group_CA2 = group_CA2[-miss,]
t_group_CA2 = t(group_CA2)
umap_results <- umap(t_group_CA2,n_neighbors=6,min_dist=0.1,n_components = 6)
umaps = data.frame(umap_results[["layout"]])
umaps$condition = as.factor(meta$condition)
ggplot(umaps,mapping = aes(x=X1,y=X2,col=condition)) +
  geom_point(aes(colour = condition))

cor_con_batch = cor(meta$condition_num,meta$batch,method = 'pearson')
cor_con_reads = cor(meta$condition_num,meta$reads,method = 'pearson')
cor_con_Align_rate = cor(meta$condition_num,meta$Align_rate,method = 'pearson')
cor_con_exon = cor(meta$condition_num,meta$exon,method = 'pearson')
cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
# -con
# exon -UMAP
corr_reads_1 = cor(umaps$X1,meta$reads,method = 'pearson')
corr_reads_2 = cor(umaps$X2,meta$reads,method = 'pearson')
corr_reads_3 = cor(umaps$X3,meta$reads,method = 'pearson')
corr_reads_4 = cor(umaps$X4,meta$reads,method = 'pearson')
corr_reads_5 = cor(umaps$X5,meta$reads,method = 'pearson')
corr_reads_6 = cor(umaps$X6,meta$reads,method = 'pearson')

corr_Align_rate_1 = cor(umaps$X1,meta$Align_rate,method = 'pearson')
corr_Align_rate_2 = cor(umaps$X2,meta$Align_rate,method = 'pearson')
corr_Align_rate_3 = cor(umaps$X3,meta$Align_rate,method = 'pearson')
corr_Align_rate_4 = cor(umaps$X4,meta$Align_rate,method = 'pearson')
corr_Align_rate_5 = cor(umaps$X5,meta$Align_rate,method = 'pearson')
corr_Align_rate_6 = cor(umaps$X6,meta$Align_rate,method = 'pearson')

corr_exon_1 = cor(umaps$X1,meta$exon,method = 'pearson')
corr_exon_2 = cor(umaps$X2,meta$exon,method = 'pearson')
corr_exon_3 = cor(umaps$X3,meta$exon,method = 'pearson')
corr_exon_4 = cor(umaps$X4,meta$exon,method = 'pearson')
corr_exon_5 = cor(umaps$X5,meta$exon,method = 'pearson')
corr_exon_6 = cor(umaps$X6,meta$exon,method = 'pearson')

corr_sequencing_depth_1 = cor(umaps$X1,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_2 = cor(umaps$X2,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_3 = cor(umaps$X3,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_4 = cor(umaps$X4,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_5 = cor(umaps$X5,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_6 = cor(umaps$X6,meta$sequencing_depth,method = 'pearson')

corr_batch_1 = cor(umaps$X1,meta$batch,method = 'pearson')
corr_batch_2 = cor(umaps$X2,meta$batch,method = 'pearson')
corr_batch_3 = cor(umaps$X3,meta$batch,method = 'pearson')
corr_batch_4 = cor(umaps$X4,meta$batch,method = 'pearson')
corr_batch_5 = cor(umaps$X5,meta$batch,method = 'pearson')
corr_batch_6 = cor(umaps$X6,meta$batch,method = 'pearson')

group_list1 = meta$condition
group_list2 = meta$reads
group_list3 = meta$Align_rate
group_list4 = meta$exon
group_list5 = meta$sequencing_depth
group_list6 = meta$batch

colData=data.frame(row.names = colnames(group_CA2),
                   condition=group_list1,
                   reads=group_list2,
                   Align_rate=group_list3,
                   exon = group_list4,
                   sequencing_depth = group_list5)

colData$condition = factor(colData$condition,levels = c('C','A'))
colData$reads =as.numeric(colData$reads)
colData$Align_rate=as.numeric(colData$Align_rate)
colData$exon=as.numeric(colData$exon)
colData$sequencing_depth=as.numeric(colData$sequencing_depth)

dds=DESeqDataSetFromMatrix(countData = group_CA2,
                           colData = colData,
                           design = ~ condition + Align_rate)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","C","A"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.6,]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_CA2[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0),]
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/RNA/DEG")
all_glia = group_CA2[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_CA2.csv')
write.csv(DEG_matrix,'Differential_CA2.csv')
###############################C vs. A_3###############################################################
meta = bind_rows(factors[43:54,1:7],factors[19:30,1:7])
meta$condition_num = c(rep(1,12),rep(2,12))
group_CA3 = data.frame(GroupC,GroupA_3)
miss_1 <- c()
for (i in 1:nrow(group_CA3)) {
  if(length(which(group_CA3[i,1:12]<10)) > 0.5*ncol(group_CA3[,1:12])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(group_CA3)) {
  if(length(which(group_CA3[i,13:24]<10)) > 0.5*ncol(group_CA3[,13:24])) miss_2 <- append(miss_2,i)
}
miss = union(miss_1,miss_2)
group_CA3 = group_CA3[-miss,]
t_group_CA3 = t(group_CA3)
umap_results <- umap(t_group_CA3,n_neighbors=12,min_dist=0.1,n_components = 6)
umaps = data.frame(umap_results[["layout"]])
umaps$condition = as.factor(meta$condition)
ggplot(umaps,mapping = aes(x=X1,y=X2,col=condition)) +
  geom_point(aes(colour = condition))

cor_con_batch = cor(meta$condition_num,meta$batch,method = 'pearson')
cor_con_reads = cor(meta$condition_num,meta$reads,method = 'pearson')
cor_con_Align_rate = cor(meta$condition_num,meta$Align_rate,method = 'pearson')
cor_con_exon = cor(meta$condition_num,meta$exon,method = 'pearson')
cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
# exon reads -con
# exon reads sequencing_depth -UMAP
corr_reads_1 = cor(umaps$X1,meta$reads,method = 'pearson')
corr_reads_2 = cor(umaps$X2,meta$reads,method = 'pearson')
corr_reads_3 = cor(umaps$X3,meta$reads,method = 'pearson')
corr_reads_4 = cor(umaps$X4,meta$reads,method = 'pearson')
corr_reads_5 = cor(umaps$X5,meta$reads,method = 'pearson')
corr_reads_6 = cor(umaps$X6,meta$reads,method = 'pearson')

corr_Align_rate_1 = cor(umaps$X1,meta$Align_rate,method = 'pearson')
corr_Align_rate_2 = cor(umaps$X2,meta$Align_rate,method = 'pearson')
corr_Align_rate_3 = cor(umaps$X3,meta$Align_rate,method = 'pearson')
corr_Align_rate_4 = cor(umaps$X4,meta$Align_rate,method = 'pearson')
corr_Align_rate_5 = cor(umaps$X5,meta$Align_rate,method = 'pearson')
corr_Align_rate_6 = cor(umaps$X6,meta$Align_rate,method = 'pearson')

corr_exon_1 = cor(umaps$X1,meta$exon,method = 'pearson')
corr_exon_2 = cor(umaps$X2,meta$exon,method = 'pearson')
corr_exon_3 = cor(umaps$X3,meta$exon,method = 'pearson')
corr_exon_4 = cor(umaps$X4,meta$exon,method = 'pearson')
corr_exon_5 = cor(umaps$X5,meta$exon,method = 'pearson')
corr_exon_6 = cor(umaps$X6,meta$exon,method = 'pearson')

corr_sequencing_depth_1 = cor(umaps$X1,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_2 = cor(umaps$X2,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_3 = cor(umaps$X3,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_4 = cor(umaps$X4,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_5 = cor(umaps$X5,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_6 = cor(umaps$X6,meta$sequencing_depth,method = 'pearson')

corr_batch_1 = cor(umaps$X1,meta$batch,method = 'pearson')
corr_batch_2 = cor(umaps$X2,meta$batch,method = 'pearson')
corr_batch_3 = cor(umaps$X3,meta$batch,method = 'pearson')
corr_batch_4 = cor(umaps$X4,meta$batch,method = 'pearson')
corr_batch_5 = cor(umaps$X5,meta$batch,method = 'pearson')
corr_batch_6 = cor(umaps$X6,meta$batch,method = 'pearson')

group_list1 = meta$condition
group_list2 = meta$reads
group_list3 = meta$Align_rate
group_list4 = meta$exon
group_list5 = meta$sequencing_depth
group_list6 = meta$batch

colData=data.frame(row.names = colnames(group_CA3),
                   condition=group_list1,
                   reads=group_list2,
                   Align_rate=group_list3,
                   exon = group_list4,
                   sequencing_depth = group_list5)

colData$condition = factor(colData$condition,levels = c('C','A'))
colData$reads =as.numeric(colData$reads)
colData$Align_rate=as.numeric(colData$Align_rate)
colData$exon=as.numeric(colData$exon)
colData$sequencing_depth=as.numeric(colData$sequencing_depth)

dds=DESeqDataSetFromMatrix(countData = group_CA3,
                           colData = colData,
                           design = ~ condition + exon + reads)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","C","A"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.6,]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_CA3[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0),]
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/RNA/DEG")
all_glia = group_CA3[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_CA3.csv')
write.csv(DEG_matrix,'Differential_CA3.csv')
###############################B vs. D_1###############################################################
meta = bind_rows(factors[31:42,1:7],factors[55:66,1:7])
meta$condition_num = c(rep(1,12),rep(2,12))
group_BD1 = data.frame(GroupB,GroupD_1)
miss_1 <- c()
for (i in 1:nrow(group_BD1)) {
  if(length(which(group_BD1[i,1:12]<10)) > 0.5*ncol(group_BD1[,1:12])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(group_BD1)) {
  if(length(which(group_BD1[i,13:24]<10)) > 0.5*ncol(group_BD1[,13:24])) miss_2 <- append(miss_2,i)
}
miss = union(miss_1,miss_2)
group_BD1 = group_BD1[-miss,]
t_group_BD1 = t(group_BD1)
umap_results <- umap(t_group_BD1,n_neighbors=12,min_dist=0.1,n_components = 6)
umaps = data.frame(umap_results[["layout"]])
umaps$condition = as.factor(meta$condition)
ggplot(umaps,mapping = aes(x=X1,y=X2,col=condition)) +
  geom_point(aes(colour = condition))

cor_con_batch = cor(meta$condition_num,meta$batch,method = 'pearson')
cor_con_reads = cor(meta$condition_num,meta$reads,method = 'pearson')
cor_con_Align_rate = cor(meta$condition_num,meta$Align_rate,method = 'pearson')
cor_con_exon = cor(meta$condition_num,meta$exon,method = 'pearson')
cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
#  -con
# reads sequencing_depth -UMAP
corr_reads_1 = cor(umaps$X1,meta$reads,method = 'pearson')
corr_reads_2 = cor(umaps$X2,meta$reads,method = 'pearson')
corr_reads_3 = cor(umaps$X3,meta$reads,method = 'pearson')
corr_reads_4 = cor(umaps$X4,meta$reads,method = 'pearson')
corr_reads_5 = cor(umaps$X5,meta$reads,method = 'pearson')
corr_reads_6 = cor(umaps$X6,meta$reads,method = 'pearson')

corr_Align_rate_1 = cor(umaps$X1,meta$Align_rate,method = 'pearson')
corr_Align_rate_2 = cor(umaps$X2,meta$Align_rate,method = 'pearson')
corr_Align_rate_3 = cor(umaps$X3,meta$Align_rate,method = 'pearson')
corr_Align_rate_4 = cor(umaps$X4,meta$Align_rate,method = 'pearson')
corr_Align_rate_5 = cor(umaps$X5,meta$Align_rate,method = 'pearson')
corr_Align_rate_6 = cor(umaps$X6,meta$Align_rate,method = 'pearson')

corr_exon_1 = cor(umaps$X1,meta$exon,method = 'pearson')
corr_exon_2 = cor(umaps$X2,meta$exon,method = 'pearson')
corr_exon_3 = cor(umaps$X3,meta$exon,method = 'pearson')
corr_exon_4 = cor(umaps$X4,meta$exon,method = 'pearson')
corr_exon_5 = cor(umaps$X5,meta$exon,method = 'pearson')
corr_exon_6 = cor(umaps$X6,meta$exon,method = 'pearson')

corr_sequencing_depth_1 = cor(umaps$X1,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_2 = cor(umaps$X2,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_3 = cor(umaps$X3,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_4 = cor(umaps$X4,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_5 = cor(umaps$X5,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_6 = cor(umaps$X6,meta$sequencing_depth,method = 'pearson')

corr_batch_1 = cor(umaps$X1,meta$batch,method = 'pearson')
corr_batch_2 = cor(umaps$X2,meta$batch,method = 'pearson')
corr_batch_3 = cor(umaps$X3,meta$batch,method = 'pearson')
corr_batch_4 = cor(umaps$X4,meta$batch,method = 'pearson')
corr_batch_5 = cor(umaps$X5,meta$batch,method = 'pearson')
corr_batch_6 = cor(umaps$X6,meta$batch,method = 'pearson')

group_list1 = meta$condition
group_list2 = meta$reads
group_list3 = meta$Align_rate
group_list4 = meta$exon
group_list5 = meta$sequencing_depth
group_list6 = meta$batch

colData=data.frame(row.names = colnames(group_BD1),
                   condition=group_list1,
                   reads=group_list2,
                   Align_rate=group_list3,
                   exon = group_list4,
                   sequencing_depth = group_list5)

colData$condition = factor(colData$condition,levels = c('B','D'))
colData$reads =as.numeric(colData$reads)
colData$Align_rate=as.numeric(colData$Align_rate)
colData$exon=as.numeric(colData$exon)
colData$sequencing_depth=as.numeric(colData$sequencing_depth)

dds=DESeqDataSetFromMatrix(countData = group_BD1,
                           colData = colData,
                           design = ~ condition + Align_rate + reads +sequencing_depth)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","B","D"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.6,]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_BD1[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0),]
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/RNA/DEG")
all_glia = group_BD1[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_BD1.csv')
write.csv(DEG_matrix,'Differential_BD1.csv')
###############################B vs. D_2###############################################################
meta = bind_rows(factors[31:42,1:7],factors[67:72,1:7])
meta$condition_num = c(rep(1,12),rep(2,6))
group_BD2 = data.frame(GroupB,GroupD_2)
miss_1 <- c()
for (i in 1:nrow(group_BD2)) {
  if(length(which(group_BD2[i,1:12]<10)) > 0.5*ncol(group_BD2[,1:12])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(group_BD2)) {
  if(length(which(group_BD2[i,13:18]<10)) > 0.5*ncol(group_BD2[,13:18])) miss_2 <- append(miss_2,i)
}
miss = union(miss_1,miss_2)
group_BD2 = group_BD2[-miss,]
t_group_BD2 = t(group_BD2)
umap_results <- umap(t_group_BD2,n_neighbors=6,min_dist=0.1,n_components = 6)
umaps = data.frame(umap_results[["layout"]])
umaps$condition = as.factor(meta$condition)
ggplot(umaps,mapping = aes(x=X1,y=X2,col=condition)) +
  geom_point(aes(colour = condition))

cor_con_batch = cor(meta$condition_num,meta$batch,method = 'pearson')
cor_con_reads = cor(meta$condition_num,meta$reads,method = 'pearson')
cor_con_Align_rate = cor(meta$condition_num,meta$Align_rate,method = 'pearson')
cor_con_exon = cor(meta$condition_num,meta$exon,method = 'pearson')
cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
#  -con
# reads sequencing_depth -UMAP
corr_reads_1 = cor(umaps$X1,meta$reads,method = 'pearson')
corr_reads_2 = cor(umaps$X2,meta$reads,method = 'pearson')
corr_reads_3 = cor(umaps$X3,meta$reads,method = 'pearson')
corr_reads_4 = cor(umaps$X4,meta$reads,method = 'pearson')
corr_reads_5 = cor(umaps$X5,meta$reads,method = 'pearson')
corr_reads_6 = cor(umaps$X6,meta$reads,method = 'pearson')

corr_Align_rate_1 = cor(umaps$X1,meta$Align_rate,method = 'pearson')
corr_Align_rate_2 = cor(umaps$X2,meta$Align_rate,method = 'pearson')
corr_Align_rate_3 = cor(umaps$X3,meta$Align_rate,method = 'pearson')
corr_Align_rate_4 = cor(umaps$X4,meta$Align_rate,method = 'pearson')
corr_Align_rate_5 = cor(umaps$X5,meta$Align_rate,method = 'pearson')
corr_Align_rate_6 = cor(umaps$X6,meta$Align_rate,method = 'pearson')

corr_exon_1 = cor(umaps$X1,meta$exon,method = 'pearson')
corr_exon_2 = cor(umaps$X2,meta$exon,method = 'pearson')
corr_exon_3 = cor(umaps$X3,meta$exon,method = 'pearson')
corr_exon_4 = cor(umaps$X4,meta$exon,method = 'pearson')
corr_exon_5 = cor(umaps$X5,meta$exon,method = 'pearson')
corr_exon_6 = cor(umaps$X6,meta$exon,method = 'pearson')

corr_sequencing_depth_1 = cor(umaps$X1,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_2 = cor(umaps$X2,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_3 = cor(umaps$X3,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_4 = cor(umaps$X4,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_5 = cor(umaps$X5,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_6 = cor(umaps$X6,meta$sequencing_depth,method = 'pearson')

corr_batch_1 = cor(umaps$X1,meta$batch,method = 'pearson')
corr_batch_2 = cor(umaps$X2,meta$batch,method = 'pearson')
corr_batch_3 = cor(umaps$X3,meta$batch,method = 'pearson')
corr_batch_4 = cor(umaps$X4,meta$batch,method = 'pearson')
corr_batch_5 = cor(umaps$X5,meta$batch,method = 'pearson')
corr_batch_6 = cor(umaps$X6,meta$batch,method = 'pearson')

group_list1 = meta$condition
group_list2 = meta$reads
group_list3 = meta$Align_rate
group_list4 = meta$exon
group_list5 = meta$sequencing_depth
group_list6 = meta$batch

colData=data.frame(row.names = colnames(group_BD2),
                   condition=group_list1,
                   reads=group_list2,
                   Align_rate=group_list3,
                   exon = group_list4,
                   sequencing_depth = group_list5)

colData$condition = factor(colData$condition,levels = c('B','D'))
colData$reads =as.numeric(colData$reads)
colData$Align_rate=as.numeric(colData$Align_rate)
colData$exon=as.numeric(colData$exon)
colData$sequencing_depth=as.numeric(colData$sequencing_depth)

dds=DESeqDataSetFromMatrix(countData = group_BD2,
                           colData = colData,
                           design = ~ condition + exon +reads)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","B","D"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.6,]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_BD2[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0),]
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/RNA/DEG")
all_glia = group_BD2[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_BD2.csv')
write.csv(DEG_matrix,'Differential_BD2.csv')
###############################B vs. D_3###############################################################
meta = bind_rows(factors[31:42,1:7],factors[73:82,1:7])
meta$condition_num = c(rep(1,12),rep(2,10))
group_BD3 = data.frame(GroupB,GroupD_3)
miss_1 <- c()
for (i in 1:nrow(group_BD3)) {
  if(length(which(group_BD3[i,1:12]<10)) > 0.5*ncol(group_BD3[,1:12])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(group_BD3)) {
  if(length(which(group_BD3[i,13:22]<10)) > 0.5*ncol(group_BD3[,13:22])) miss_2 <- append(miss_2,i)
}
miss = union(miss_1,miss_2)
group_BD3 = group_BD3[-miss,]
t_group_BD3 = t(group_BD3)
umap_results <- umap(t_group_BD3,n_neighbors=10,min_dist=0.1,n_components = 6)
umaps = data.frame(umap_results[["layout"]])
umaps$condition = as.factor(meta$condition)
ggplot(umaps,mapping = aes(x=X1,y=X2,col=condition)) +
  geom_point(aes(colour = condition))

cor_con_batch = cor(meta$condition_num,meta$batch,method = 'pearson')
cor_con_reads = cor(meta$condition_num,meta$reads,method = 'pearson')
cor_con_Align_rate = cor(meta$condition_num,meta$Align_rate,method = 'pearson')
cor_con_exon = cor(meta$condition_num,meta$exon,method = 'pearson')
cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
#  -con
# reads sequencing_depth -UMAP
corr_reads_1 = cor(umaps$X1,meta$reads,method = 'pearson')
corr_reads_2 = cor(umaps$X2,meta$reads,method = 'pearson')
corr_reads_3 = cor(umaps$X3,meta$reads,method = 'pearson')
corr_reads_4 = cor(umaps$X4,meta$reads,method = 'pearson')
corr_reads_5 = cor(umaps$X5,meta$reads,method = 'pearson')
corr_reads_6 = cor(umaps$X6,meta$reads,method = 'pearson')

corr_Align_rate_1 = cor(umaps$X1,meta$Align_rate,method = 'pearson')
corr_Align_rate_2 = cor(umaps$X2,meta$Align_rate,method = 'pearson')
corr_Align_rate_3 = cor(umaps$X3,meta$Align_rate,method = 'pearson')
corr_Align_rate_4 = cor(umaps$X4,meta$Align_rate,method = 'pearson')
corr_Align_rate_5 = cor(umaps$X5,meta$Align_rate,method = 'pearson')
corr_Align_rate_6 = cor(umaps$X6,meta$Align_rate,method = 'pearson')

corr_exon_1 = cor(umaps$X1,meta$exon,method = 'pearson')
corr_exon_2 = cor(umaps$X2,meta$exon,method = 'pearson')
corr_exon_3 = cor(umaps$X3,meta$exon,method = 'pearson')
corr_exon_4 = cor(umaps$X4,meta$exon,method = 'pearson')
corr_exon_5 = cor(umaps$X5,meta$exon,method = 'pearson')
corr_exon_6 = cor(umaps$X6,meta$exon,method = 'pearson')

corr_sequencing_depth_1 = cor(umaps$X1,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_2 = cor(umaps$X2,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_3 = cor(umaps$X3,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_4 = cor(umaps$X4,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_5 = cor(umaps$X5,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_6 = cor(umaps$X6,meta$sequencing_depth,method = 'pearson')

corr_batch_1 = cor(umaps$X1,meta$batch,method = 'pearson')
corr_batch_2 = cor(umaps$X2,meta$batch,method = 'pearson')
corr_batch_3 = cor(umaps$X3,meta$batch,method = 'pearson')
corr_batch_4 = cor(umaps$X4,meta$batch,method = 'pearson')
corr_batch_5 = cor(umaps$X5,meta$batch,method = 'pearson')
corr_batch_6 = cor(umaps$X6,meta$batch,method = 'pearson')

group_list1 = meta$condition
group_list2 = meta$reads
group_list3 = meta$Align_rate
group_list4 = meta$exon
group_list5 = meta$sequencing_depth
group_list6 = meta$batch

colData=data.frame(row.names = colnames(group_BD3),
                   condition=group_list1,
                   reads=group_list2,
                   Align_rate=group_list3,
                   exon = group_list4,
                   sequencing_depth = group_list5)

colData$condition = factor(colData$condition,levels = c('B','D'))
colData$reads =as.numeric(colData$reads)
colData$Align_rate=as.numeric(colData$Align_rate)
colData$exon=as.numeric(colData$exon)
colData$sequencing_depth=as.numeric(colData$sequencing_depth)

dds=DESeqDataSetFromMatrix(countData = group_BD3,
                           colData = colData,
                           design = ~ condition + exon + reads)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","B","D"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.6,]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_BD3[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0),]
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/RNA/DEG")
all_glia = group_BD3[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_BD3.csv')
write.csv(DEG_matrix,'Differential_BD3.csv')
###############################F vs. G###############################################################
meta = bind_rows(factors[93:102,1:7],factors[103:114,1:7])
meta$condition_num = c(rep(1,10),rep(2,12))
group_FG = data.frame(GroupF,GroupG)
miss_1 <- c()
for (i in 1:nrow(group_FG)) {
  if(length(which(group_FG[i,1:10]<10)) > 0.5*ncol(group_FG[,1:10])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(group_FG)) {
  if(length(which(group_FG[i,11:22]<10)) > 0.5*ncol(group_FG[,11:22])) miss_2 <- append(miss_2,i)
}
miss = union(miss_1,miss_2)
group_FG = group_FG[-miss,]
t_group_FG = t(group_FG)
umap_results <- umap(t_group_FG,n_neighbors=10,min_dist=0.1,n_components = 6)
umaps = data.frame(umap_results[["layout"]])
umaps$condition = as.factor(meta$condition)
ggplot(umaps,mapping = aes(x=X1,y=X2,col=condition)) +
  geom_point(aes(colour = condition))

cor_con_batch = cor(meta$condition_num,meta$batch,method = 'pearson')
cor_con_reads = cor(meta$condition_num,meta$reads,method = 'pearson')
cor_con_Align_rate = cor(meta$condition_num,meta$Align_rate,method = 'pearson')
cor_con_exon = cor(meta$condition_num,meta$exon,method = 'pearson')
cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
# Aligned_rate exon  -con
# Aligned_rate exon reads sequencing_depth -UMAP
corr_reads_1 = cor(umaps$X1,meta$reads,method = 'pearson')
corr_reads_2 = cor(umaps$X2,meta$reads,method = 'pearson')
corr_reads_3 = cor(umaps$X3,meta$reads,method = 'pearson')
corr_reads_4 = cor(umaps$X4,meta$reads,method = 'pearson')
corr_reads_5 = cor(umaps$X5,meta$reads,method = 'pearson')
corr_reads_6 = cor(umaps$X6,meta$reads,method = 'pearson')

corr_Align_rate_1 = cor(umaps$X1,meta$Align_rate,method = 'pearson')
corr_Align_rate_2 = cor(umaps$X2,meta$Align_rate,method = 'pearson')
corr_Align_rate_3 = cor(umaps$X3,meta$Align_rate,method = 'pearson')
corr_Align_rate_4 = cor(umaps$X4,meta$Align_rate,method = 'pearson')
corr_Align_rate_5 = cor(umaps$X5,meta$Align_rate,method = 'pearson')
corr_Align_rate_6 = cor(umaps$X6,meta$Align_rate,method = 'pearson')

corr_exon_1 = cor(umaps$X1,meta$exon,method = 'pearson')
corr_exon_2 = cor(umaps$X2,meta$exon,method = 'pearson')
corr_exon_3 = cor(umaps$X3,meta$exon,method = 'pearson')
corr_exon_4 = cor(umaps$X4,meta$exon,method = 'pearson')
corr_exon_5 = cor(umaps$X5,meta$exon,method = 'pearson')
corr_exon_6 = cor(umaps$X6,meta$exon,method = 'pearson')

corr_sequencing_depth_1 = cor(umaps$X1,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_2 = cor(umaps$X2,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_3 = cor(umaps$X3,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_4 = cor(umaps$X4,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_5 = cor(umaps$X5,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_6 = cor(umaps$X6,meta$sequencing_depth,method = 'pearson')

corr_batch_1 = cor(umaps$X1,meta$batch,method = 'pearson')
corr_batch_2 = cor(umaps$X2,meta$batch,method = 'pearson')
corr_batch_3 = cor(umaps$X3,meta$batch,method = 'pearson')
corr_batch_4 = cor(umaps$X4,meta$batch,method = 'pearson')
corr_batch_5 = cor(umaps$X5,meta$batch,method = 'pearson')
corr_batch_6 = cor(umaps$X6,meta$batch,method = 'pearson')

group_list1 = meta$condition
group_list2 = meta$reads
group_list3 = meta$Align_rate
group_list4 = meta$exon
group_list5 = meta$sequencing_depth
group_list6 = meta$batch

colData=data.frame(row.names = colnames(group_FG),
                   condition=group_list1,
                   reads=group_list2,
                   Align_rate=group_list3,
                   exon = group_list4,
                   sequencing_depth = group_list5)

colData$condition = factor(colData$condition,levels = c('F','G'))
colData$reads =as.numeric(colData$reads)
colData$Align_rate=as.numeric(colData$Align_rate)
colData$exon=as.numeric(colData$exon)
colData$sequencing_depth=as.numeric(colData$sequencing_depth)

dds=DESeqDataSetFromMatrix(countData = group_FG,
                           colData = colData,
                           design = ~ condition + Align_rate + exon)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","F","G"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.6,]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_FG[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0),]
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/RNA/DEG")
all_glia = group_FG[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_FG.csv')
write.csv(DEG_matrix,'Differential_FG.csv')

###############################E vs. G###############################################################
meta = bind_rows(factors[83:92,1:7],factors[103:114,1:7])
meta$condition_num = c(rep(1,10),rep(2,12))
group_EG = data.frame(GroupE,GroupG)
miss_1 <- c()
for (i in 1:nrow(group_EG)) {
  if(length(which(group_EG[i,1:12]<10)) > 0.5*ncol(group_EG[,1:12])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(group_EG)) {
  if(length(which(group_EG[i,13:22]<10)) > 0.5*ncol(group_EG[,13:22])) miss_2 <- append(miss_2,i)
}
miss = union(miss_1,miss_2)
group_EG = group_EG[-miss,]
t_group_EG = t(group_EG)
umap_results <- umap(t_group_EG,n_neighbors=10,min_dist=0.1,n_components = 6)
umaps = data.frame(umap_results[["layout"]])
umaps$condition = as.factor(meta$condition)
ggplot(umaps,mapping = aes(x=X1,y=X2,col=condition)) +
  geom_point(aes(colour = condition))

cor_con_batch = cor(meta$condition_num,meta$batch,method = 'pearson')
cor_con_reads = cor(meta$condition_num,meta$reads,method = 'pearson')
cor_con_Align_rate = cor(meta$condition_num,meta$Align_rate,method = 'pearson')
cor_con_exon = cor(meta$condition_num,meta$exon,method = 'pearson')
cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
# -con
# exon reads sequencing_depth -UMAP
corr_reads_1 = cor(umaps$X1,meta$reads,method = 'pearson')
corr_reads_2 = cor(umaps$X2,meta$reads,method = 'pearson')
corr_reads_3 = cor(umaps$X3,meta$reads,method = 'pearson')
corr_reads_4 = cor(umaps$X4,meta$reads,method = 'pearson')
corr_reads_5 = cor(umaps$X5,meta$reads,method = 'pearson')
corr_reads_6 = cor(umaps$X6,meta$reads,method = 'pearson')

corr_Align_rate_1 = cor(umaps$X1,meta$Align_rate,method = 'pearson')
corr_Align_rate_2 = cor(umaps$X2,meta$Align_rate,method = 'pearson')
corr_Align_rate_3 = cor(umaps$X3,meta$Align_rate,method = 'pearson')
corr_Align_rate_4 = cor(umaps$X4,meta$Align_rate,method = 'pearson')
corr_Align_rate_5 = cor(umaps$X5,meta$Align_rate,method = 'pearson')
corr_Align_rate_6 = cor(umaps$X6,meta$Align_rate,method = 'pearson')

corr_exon_1 = cor(umaps$X1,meta$exon,method = 'pearson')
corr_exon_2 = cor(umaps$X2,meta$exon,method = 'pearson')
corr_exon_3 = cor(umaps$X3,meta$exon,method = 'pearson')
corr_exon_4 = cor(umaps$X4,meta$exon,method = 'pearson')
corr_exon_5 = cor(umaps$X5,meta$exon,method = 'pearson')
corr_exon_6 = cor(umaps$X6,meta$exon,method = 'pearson')

corr_sequencing_depth_1 = cor(umaps$X1,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_2 = cor(umaps$X2,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_3 = cor(umaps$X3,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_4 = cor(umaps$X4,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_5 = cor(umaps$X5,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_6 = cor(umaps$X6,meta$sequencing_depth,method = 'pearson')

corr_batch_1 = cor(umaps$X1,meta$batch,method = 'pearson')
corr_batch_2 = cor(umaps$X2,meta$batch,method = 'pearson')
corr_batch_3 = cor(umaps$X3,meta$batch,method = 'pearson')
corr_batch_4 = cor(umaps$X4,meta$batch,method = 'pearson')
corr_batch_5 = cor(umaps$X5,meta$batch,method = 'pearson')
corr_batch_6 = cor(umaps$X6,meta$batch,method = 'pearson')

group_list1 = meta$condition
group_list2 = meta$reads
group_list3 = meta$Align_rate
group_list4 = meta$exon
group_list5 = meta$sequencing_depth
group_list6 = meta$batch

colData=data.frame(row.names = colnames(group_EG),
                   condition=group_list1,
                   reads=group_list2,
                   Align_rate=group_list3,
                   exon = group_list4,
                   sequencing_depth = group_list5)

colData$condition = factor(colData$condition,levels = c('E','G'))
colData$reads =as.numeric(colData$reads)
colData$Align_rate=as.numeric(colData$Align_rate)
colData$exon=as.numeric(colData$exon)
colData$sequencing_depth=as.numeric(colData$sequencing_depth)

dds=DESeqDataSetFromMatrix(countData = group_EG,
                           colData = colData,
                           design = ~ condition)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","E","G"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.6,]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_EG[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0),]
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/RNA/DEG")
all_glia = group_EG[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_EG.csv')
write.csv(DEG_matrix,'Differential_EG.csv')
###############################F vs. E###############################################################
meta = bind_rows(factors[93:102,1:7],factors[83:92,1:7])
meta$condition_num = c(rep(1,10),rep(2,10))
group_FE = data.frame(GroupF,GroupE)
miss_1 <- c()
for (i in 1:nrow(group_FE)) {
  if(length(which(group_FE[i,1:10]<10)) > 0.5*ncol(group_FE[,1:10])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(group_FE)) {
  if(length(which(group_FE[i,11:20]<10)) > 0.5*ncol(group_FE[,11:20])) miss_2 <- append(miss_2,i)
}
miss = union(miss_1,miss_2)
group_FE = group_FE[-miss,]
t_group_FE = t(group_FE)
umap_results <- umap(t_group_FE,n_neighbors=10,min_dist=0.1,n_components = 6)
umaps = data.frame(umap_results[["layout"]])
umaps$condition = as.factor(meta$condition)
ggplot(umaps,mapping = aes(x=X1,y=X2,col=condition)) +
  geom_point(aes(colour = condition))

cor_con_batch = cor(meta$condition_num,meta$batch,method = 'pearson')
cor_con_reads = cor(meta$condition_num,meta$reads,method = 'pearson')
cor_con_Align_rate = cor(meta$condition_num,meta$Align_rate,method = 'pearson')
cor_con_exon = cor(meta$condition_num,meta$exon,method = 'pearson')
cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
# exon reads -con
# Aligned_rate exon reads sequencing_depth -UMAP
corr_reads_1 = cor(umaps$X1,meta$reads,method = 'pearson')
corr_reads_2 = cor(umaps$X2,meta$reads,method = 'pearson')
corr_reads_3 = cor(umaps$X3,meta$reads,method = 'pearson')
corr_reads_4 = cor(umaps$X4,meta$reads,method = 'pearson')
corr_reads_5 = cor(umaps$X5,meta$reads,method = 'pearson')
corr_reads_6 = cor(umaps$X6,meta$reads,method = 'pearson')

corr_Align_rate_1 = cor(umaps$X1,meta$Align_rate,method = 'pearson')
corr_Align_rate_2 = cor(umaps$X2,meta$Align_rate,method = 'pearson')
corr_Align_rate_3 = cor(umaps$X3,meta$Align_rate,method = 'pearson')
corr_Align_rate_4 = cor(umaps$X4,meta$Align_rate,method = 'pearson')
corr_Align_rate_5 = cor(umaps$X5,meta$Align_rate,method = 'pearson')
corr_Align_rate_6 = cor(umaps$X6,meta$Align_rate,method = 'pearson')

corr_exon_1 = cor(umaps$X1,meta$exon,method = 'pearson')
corr_exon_2 = cor(umaps$X2,meta$exon,method = 'pearson')
corr_exon_3 = cor(umaps$X3,meta$exon,method = 'pearson')
corr_exon_4 = cor(umaps$X4,meta$exon,method = 'pearson')
corr_exon_5 = cor(umaps$X5,meta$exon,method = 'pearson')
corr_exon_6 = cor(umaps$X6,meta$exon,method = 'pearson')

corr_sequencing_depth_1 = cor(umaps$X1,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_2 = cor(umaps$X2,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_3 = cor(umaps$X3,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_4 = cor(umaps$X4,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_5 = cor(umaps$X5,meta$sequencing_depth,method = 'pearson')
corr_sequencing_depth_6 = cor(umaps$X6,meta$sequencing_depth,method = 'pearson')

corr_batch_1 = cor(umaps$X1,meta$batch,method = 'pearson')
corr_batch_2 = cor(umaps$X2,meta$batch,method = 'pearson')
corr_batch_3 = cor(umaps$X3,meta$batch,method = 'pearson')
corr_batch_4 = cor(umaps$X4,meta$batch,method = 'pearson')
corr_batch_5 = cor(umaps$X5,meta$batch,method = 'pearson')
corr_batch_6 = cor(umaps$X6,meta$batch,method = 'pearson')

group_list1 = meta$condition
group_list2 = meta$reads
group_list3 = meta$Align_rate
group_list4 = meta$exon
group_list5 = meta$sequencing_depth
group_list6 = meta$batch

colData=data.frame(row.names = colnames(group_FE),
                   condition=group_list1,
                   reads=group_list2,
                   Align_rate=group_list3,
                   exon = group_list4,
                   sequencing_depth = group_list5,
                   batch = group_list6)

colData$condition = factor(colData$condition,levels = c('F','E'))
colData$reads =as.numeric(colData$reads)
colData$Align_rate=as.numeric(colData$Align_rate)
colData$exon=as.numeric(colData$exon)
colData$sequencing_depth=as.numeric(colData$sequencing_depth)
colData$batch = as.factor(colData$batch)

dds=DESeqDataSetFromMatrix(countData = group_FE,
                           colData = colData,
                           design = ~ condition + Align_rate)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","F","E"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.6,]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_FE[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0),]
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/RNA/DEG")
all_glia = group_FE[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_FE.csv')
write.csv(DEG_matrix,'Differential_FE.csv')

# ###############################D_1 vs. A_1###############################################################
# meta = bind_rows(factors[55:66,1:7],factors[1:12,1:7])
# meta$condition_num = c(rep(1,12),rep(2,12))
# group_D1A1 = data.frame(GroupD_1,GroupA_1)
# miss_1 <- c()
# for (i in 1:nrow(group_D1A1)) {
#   if(length(which(group_D1A1[i,1:12]<10)) > 0.5*ncol(group_D1A1[,1:12])) miss_1 <- append(miss_1,i)
# }
# miss_2 <- c()
# for (i in 1:nrow(group_D1A1)) {
#   if(length(which(group_D1A1[i,13:24]<10)) > 0.5*ncol(group_D1A1[,13:24])) miss_2 <- append(miss_2,i)
# }
# miss = union(miss_1,miss_2)
# group_D1A1 = group_D1A1[-miss,]
# t_group_D1A1 = t(group_D1A1)
# umap_results <- umap(t_group_D1A1,n_neighbors=12,min_dist=0.1,n_components = 6)
# umaps = data.frame(umap_results[["layout"]])
# umaps$condition = as.factor(meta$condition)
# ggplot(umaps,mapping = aes(x=X1,y=X2,col=condition)) +
#   geom_point(aes(colour = condition))
# 
# cor_con_batch = cor(meta$condition_num,meta$batch,method = 'pearson')
# cor_con_reads = cor(meta$condition_num,meta$reads,method = 'pearson')
# cor_con_Align_rate = cor(meta$condition_num,meta$Align_rate,method = 'pearson')
# cor_con_exon = cor(meta$condition_num,meta$exon,method = 'pearson')
# cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
# #  -con
# # Align_rate exon reads sequencing_depth -UMAP
# corr_reads_1 = cor(umaps$X1,meta$reads,method = 'pearson')
# corr_reads_2 = cor(umaps$X2,meta$reads,method = 'pearson')
# corr_reads_3 = cor(umaps$X3,meta$reads,method = 'pearson')
# corr_reads_4 = cor(umaps$X4,meta$reads,method = 'pearson')
# corr_reads_5 = cor(umaps$X5,meta$reads,method = 'pearson')
# corr_reads_6 = cor(umaps$X6,meta$reads,method = 'pearson')
# 
# corr_Align_rate_1 = cor(umaps$X1,meta$Align_rate,method = 'pearson')
# corr_Align_rate_2 = cor(umaps$X2,meta$Align_rate,method = 'pearson')
# corr_Align_rate_3 = cor(umaps$X3,meta$Align_rate,method = 'pearson')
# corr_Align_rate_4 = cor(umaps$X4,meta$Align_rate,method = 'pearson')
# corr_Align_rate_5 = cor(umaps$X5,meta$Align_rate,method = 'pearson')
# corr_Align_rate_6 = cor(umaps$X6,meta$Align_rate,method = 'pearson')
# 
# corr_exon_1 = cor(umaps$X1,meta$exon,method = 'pearson')
# corr_exon_2 = cor(umaps$X2,meta$exon,method = 'pearson')
# corr_exon_3 = cor(umaps$X3,meta$exon,method = 'pearson')
# corr_exon_4 = cor(umaps$X4,meta$exon,method = 'pearson')
# corr_exon_5 = cor(umaps$X5,meta$exon,method = 'pearson')
# corr_exon_6 = cor(umaps$X6,meta$exon,method = 'pearson')
# 
# corr_sequencing_depth_1 = cor(umaps$X1,meta$sequencing_depth,method = 'pearson')
# corr_sequencing_depth_2 = cor(umaps$X2,meta$sequencing_depth,method = 'pearson')
# corr_sequencing_depth_3 = cor(umaps$X3,meta$sequencing_depth,method = 'pearson')
# corr_sequencing_depth_4 = cor(umaps$X4,meta$sequencing_depth,method = 'pearson')
# corr_sequencing_depth_5 = cor(umaps$X5,meta$sequencing_depth,method = 'pearson')
# corr_sequencing_depth_6 = cor(umaps$X6,meta$sequencing_depth,method = 'pearson')
# 
# corr_batch_1 = cor(umaps$X1,meta$batch,method = 'pearson')
# corr_batch_2 = cor(umaps$X2,meta$batch,method = 'pearson')
# corr_batch_3 = cor(umaps$X3,meta$batch,method = 'pearson')
# corr_batch_4 = cor(umaps$X4,meta$batch,method = 'pearson')
# corr_batch_5 = cor(umaps$X5,meta$batch,method = 'pearson')
# corr_batch_6 = cor(umaps$X6,meta$batch,method = 'pearson')
# 
# group_list1 = meta$condition
# group_list2 = meta$reads
# group_list3 = meta$Align_rate
# group_list4 = meta$exon
# group_list5 = meta$sequencing_depth
# 
# colData=data.frame(row.names = colnames(group_D1A1),
#                    condition=group_list1,
#                    reads=group_list2,
#                    Align_rate=group_list3,
#                    exon = group_list4,
#                    sequencing_depth = group_list5)
# 
# colData$condition = factor(colData$condition,levels = c('D','A'))
# colData$reads =as.numeric(colData$reads)
# colData$Align_rate=as.numeric(colData$Align_rate)
# colData$exon=as.numeric(colData$exon)
# colData$sequencing_depth=as.numeric(colData$sequencing_depth)
# 
# dds=DESeqDataSetFromMatrix(countData = group_D1A1,
#                            colData = colData,
#                            design = ~ condition + exon + sequencing_depth + reads)
# #tmp = model.matrix(~group_list + gender + Age + PMI, colData)
# 
# dds=DESeq(dds)
# 
# res_glia=results(dds,
#                  contrast = c("condition","D","A"))
# resOrdered_glia=res_glia[order(res_glia$padj),]
# DEG_glia=as.data.frame(resOrdered_glia)
# DEG_glia=na.omit(DEG_glia)
# nrDEG_glia=DEG_glia[DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 1,]
# choose_gene_neuron=rownames(nrDEG_glia)
# choose_matrix=group_D1A1[choose_gene_neuron,]
# DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
# up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0),]
# down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0),]
# setwd("~/Documents/Virginia Tech/Project_3_MIA mouse/differential analysis/results/DEGs/ABCD/union_FC_0.6")
# all_glia = group_D1A1[rownames(DEG_glia),]
# all_matrix = data.frame(DEG_glia,all_glia)
# write.csv(all_matrix,'All_D1A1.csv')
# write.csv(DEG_matrix,'Differential_D1A1.csv')
# 
# ###############################D_2 vs. A_2###############################################################
# meta = bind_rows(factors[67:72,1:7],factors[13:18,1:7])
# meta$condition_num = c(rep(1,6),rep(2,6))
# group_D2A2 = data.frame(GroupD_2,GroupA_2)
# miss_1 <- c()
# for (i in 1:nrow(group_D2A2)) {
#   if(length(which(group_D2A2[i,1:6]<10)) > 0.5*ncol(group_D2A2[,1:6])) miss_1 <- append(miss_1,i)
# }
# miss_2 <- c()
# for (i in 1:nrow(group_D2A2)) {
#   if(length(which(group_D2A2[i,7:12]<10)) > 0.5*ncol(group_D2A2[,7:12])) miss_2 <- append(miss_2,i)
# }
# miss = union(miss_1,miss_2)
# group_D2A2 = group_D2A2[-miss,]
# t_group_D2A2 = t(group_D2A2)
# umap_results <- umap(t_group_D2A2,n_neighbors=6,min_dist=0.1,n_components = 6)
# umaps = data.frame(umap_results[["layout"]])
# umaps$condition = as.factor(meta$condition)
# ggplot(umaps,mapping = aes(x=X1,y=X2,col=condition)) +
#   geom_point(aes(colour = condition))
# 
# cor_con_batch = cor(meta$condition_num,meta$batch,method = 'pearson')
# cor_con_reads = cor(meta$condition_num,meta$reads,method = 'pearson')
# cor_con_Align_rate = cor(meta$condition_num,meta$Align_rate,method = 'pearson')
# cor_con_exon = cor(meta$condition_num,meta$exon,method = 'pearson')
# cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
# #  -con
# # -UMAP
# corr_reads_1 = cor(umaps$X1,meta$reads,method = 'pearson')
# corr_reads_2 = cor(umaps$X2,meta$reads,method = 'pearson')
# corr_reads_3 = cor(umaps$X3,meta$reads,method = 'pearson')
# corr_reads_4 = cor(umaps$X4,meta$reads,method = 'pearson')
# corr_reads_5 = cor(umaps$X5,meta$reads,method = 'pearson')
# corr_reads_6 = cor(umaps$X6,meta$reads,method = 'pearson')
# 
# corr_Align_rate_1 = cor(umaps$X1,meta$Align_rate,method = 'pearson')
# corr_Align_rate_2 = cor(umaps$X2,meta$Align_rate,method = 'pearson')
# corr_Align_rate_3 = cor(umaps$X3,meta$Align_rate,method = 'pearson')
# corr_Align_rate_4 = cor(umaps$X4,meta$Align_rate,method = 'pearson')
# corr_Align_rate_5 = cor(umaps$X5,meta$Align_rate,method = 'pearson')
# corr_Align_rate_6 = cor(umaps$X6,meta$Align_rate,method = 'pearson')
# 
# corr_exon_1 = cor(umaps$X1,meta$exon,method = 'pearson')
# corr_exon_2 = cor(umaps$X2,meta$exon,method = 'pearson')
# corr_exon_3 = cor(umaps$X3,meta$exon,method = 'pearson')
# corr_exon_4 = cor(umaps$X4,meta$exon,method = 'pearson')
# corr_exon_5 = cor(umaps$X5,meta$exon,method = 'pearson')
# corr_exon_6 = cor(umaps$X6,meta$exon,method = 'pearson')
# 
# corr_sequencing_depth_1 = cor(umaps$X1,meta$sequencing_depth,method = 'pearson')
# corr_sequencing_depth_2 = cor(umaps$X2,meta$sequencing_depth,method = 'pearson')
# corr_sequencing_depth_3 = cor(umaps$X3,meta$sequencing_depth,method = 'pearson')
# corr_sequencing_depth_4 = cor(umaps$X4,meta$sequencing_depth,method = 'pearson')
# corr_sequencing_depth_5 = cor(umaps$X5,meta$sequencing_depth,method = 'pearson')
# corr_sequencing_depth_6 = cor(umaps$X6,meta$sequencing_depth,method = 'pearson')
# 
# corr_batch_1 = cor(umaps$X1,meta$batch,method = 'pearson')
# corr_batch_2 = cor(umaps$X2,meta$batch,method = 'pearson')
# corr_batch_3 = cor(umaps$X3,meta$batch,method = 'pearson')
# corr_batch_4 = cor(umaps$X4,meta$batch,method = 'pearson')
# corr_batch_5 = cor(umaps$X5,meta$batch,method = 'pearson')
# corr_batch_6 = cor(umaps$X6,meta$batch,method = 'pearson')
# 
# group_list1 = meta$condition
# group_list2 = meta$reads
# group_list3 = meta$Align_rate
# group_list4 = meta$exon
# group_list5 = meta$sequencing_depth
# group_list6 = meta$batch
# 
# colData=data.frame(row.names = colnames(group_D2A2),
#                    condition=group_list1,
#                    reads=group_list2,
#                    Align_rate=group_list3,
#                    exon = group_list4,
#                    sequencing_depth = group_list5,
#                    batch = group_list6)
# 
# colData$condition = factor(colData$condition,levels = c('D','A'))
# colData$reads =as.numeric(colData$reads)
# colData$Align_rate=as.numeric(colData$Align_rate)
# colData$exon=as.numeric(colData$exon)
# colData$sequencing_depth=as.numeric(colData$sequencing_depth)
# colData$batch=as.factor(colData$batch)
# 
# dds=DESeqDataSetFromMatrix(countData = group_D2A2,
#                            colData = colData,
#                            design = ~ condition)
# #tmp = model.matrix(~group_list + gender + Age + PMI, colData)
# 
# dds=DESeq(dds)
# 
# res_glia=results(dds,
#                  contrast = c("condition","D","A"))
# resOrdered_glia=res_glia[order(res_glia$padj),]
# DEG_glia=as.data.frame(resOrdered_glia)
# DEG_glia=na.omit(DEG_glia)
# nrDEG_glia=DEG_glia[DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.6,]
# choose_gene_neuron=rownames(nrDEG_glia)
# choose_matrix=group_D2A2[choose_gene_neuron,]
# DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
# up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0),]
# down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0),]
# setwd("~/Documents/Virginia Tech/Project_3_MIA mouse/differential analysis/results/DEGs/ABCD/union_FC_0.6")
# all_glia = group_D2A2[rownames(DEG_glia),]
# all_matrix = data.frame(DEG_glia,all_glia)
# write.csv(all_matrix,'All_D2A2.csv')
# write.csv(DEG_matrix,'Differential_D2A2.csv')
# ###############################D_3 vs. A_3###############################################################
# meta = bind_rows(factors[73:82,1:7],factors[19:30,1:7])
# meta$condition_num = c(rep(1,10),rep(2,12))
# group_D3A3 = data.frame(GroupD_3,GroupA_3)
# miss_1 <- c()
# for (i in 1:nrow(group_D3A3)) {
#   if(length(which(group_D3A3[i,1:10]<10)) > 0.5*ncol(group_D3A3[,1:10])) miss_1 <- append(miss_1,i)
# }
# miss_2 <- c()
# for (i in 1:nrow(group_D3A3)) {
#   if(length(which(group_D3A3[i,11:22]<10)) > 0.5*ncol(group_D3A3[,11:22])) miss_2 <- append(miss_2,i)
# }
# miss = union(miss_1,miss_2)
# group_D3A3 = group_D3A3[-miss,]
# t_group_D3A3 = t(group_D3A3)
# umap_results <- umap(t_group_D3A3,n_neighbors=10,min_dist=0.1,n_components = 6)
# umaps = data.frame(umap_results[["layout"]])
# umaps$condition = as.factor(meta$condition)
# ggplot(umaps,mapping = aes(x=X1,y=X2,col=condition)) +
#   geom_point(aes(colour = condition))
# 
# cor_con_batch = cor(meta$condition_num,meta$batch,method = 'pearson')
# cor_con_reads = cor(meta$condition_num,meta$reads,method = 'pearson')
# cor_con_Align_rate = cor(meta$condition_num,meta$Align_rate,method = 'pearson')
# cor_con_exon = cor(meta$condition_num,meta$exon,method = 'pearson')
# cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
# # -con
# # reads sequencing_depth -UMAP
# corr_reads_1 = cor(umaps$X1,meta$reads,method = 'pearson')
# corr_reads_2 = cor(umaps$X2,meta$reads,method = 'pearson')
# corr_reads_3 = cor(umaps$X3,meta$reads,method = 'pearson')
# corr_reads_4 = cor(umaps$X4,meta$reads,method = 'pearson')
# corr_reads_5 = cor(umaps$X5,meta$reads,method = 'pearson')
# corr_reads_6 = cor(umaps$X6,meta$reads,method = 'pearson')
# 
# corr_Align_rate_1 = cor(umaps$X1,meta$Align_rate,method = 'pearson')
# corr_Align_rate_2 = cor(umaps$X2,meta$Align_rate,method = 'pearson')
# corr_Align_rate_3 = cor(umaps$X3,meta$Align_rate,method = 'pearson')
# corr_Align_rate_4 = cor(umaps$X4,meta$Align_rate,method = 'pearson')
# corr_Align_rate_5 = cor(umaps$X5,meta$Align_rate,method = 'pearson')
# corr_Align_rate_6 = cor(umaps$X6,meta$Align_rate,method = 'pearson')
# 
# corr_exon_1 = cor(umaps$X1,meta$exon,method = 'pearson')
# corr_exon_2 = cor(umaps$X2,meta$exon,method = 'pearson')
# corr_exon_3 = cor(umaps$X3,meta$exon,method = 'pearson')
# corr_exon_4 = cor(umaps$X4,meta$exon,method = 'pearson')
# corr_exon_5 = cor(umaps$X5,meta$exon,method = 'pearson')
# corr_exon_6 = cor(umaps$X6,meta$exon,method = 'pearson')
# 
# corr_sequencing_depth_1 = cor(umaps$X1,meta$sequencing_depth,method = 'pearson')
# corr_sequencing_depth_2 = cor(umaps$X2,meta$sequencing_depth,method = 'pearson')
# corr_sequencing_depth_3 = cor(umaps$X3,meta$sequencing_depth,method = 'pearson')
# corr_sequencing_depth_4 = cor(umaps$X4,meta$sequencing_depth,method = 'pearson')
# corr_sequencing_depth_5 = cor(umaps$X5,meta$sequencing_depth,method = 'pearson')
# corr_sequencing_depth_6 = cor(umaps$X6,meta$sequencing_depth,method = 'pearson')
# 
# corr_batch_1 = cor(umaps$X1,meta$batch,method = 'pearson')
# corr_batch_2 = cor(umaps$X2,meta$batch,method = 'pearson')
# corr_batch_3 = cor(umaps$X3,meta$batch,method = 'pearson')
# corr_batch_4 = cor(umaps$X4,meta$batch,method = 'pearson')
# corr_batch_5 = cor(umaps$X5,meta$batch,method = 'pearson')
# corr_batch_6 = cor(umaps$X6,meta$batch,method = 'pearson')
# 
# group_list1 = meta$condition
# group_list2 = meta$reads
# group_list3 = meta$Align_rate
# group_list4 = meta$exon
# group_list5 = meta$sequencing_depth
# group_list6 = meta$batch
# 
# colData=data.frame(row.names = colnames(group_D3A3),
#                    condition=group_list1,
#                    reads=group_list2,
#                    Align_rate=group_list3,
#                    exon = group_list4,
#                    sequencing_depth = group_list5,
#                    batch = group_list6)
# 
# colData$condition = factor(colData$condition,levels = c('D','A'))
# colData$reads =as.numeric(colData$reads)
# colData$Align_rate=as.numeric(colData$Align_rate)
# colData$exon=as.numeric(colData$exon)
# colData$sequencing_depth=as.numeric(colData$sequencing_depth)
# colData$batch=as.factor(colData$batch)
# 
# dds=DESeqDataSetFromMatrix(countData = group_D3A3,
#                            colData = colData,
#                            design = ~ condition + reads + sequencing_depth)
# #tmp = model.matrix(~group_list + gender + Age + PMI, colData)
# 
# dds=DESeq(dds)
# 
# res_glia=results(dds,
#                  contrast = c("condition","D","A"))
# resOrdered_glia=res_glia[order(res_glia$padj),]
# DEG_glia=as.data.frame(resOrdered_glia)
# DEG_glia=na.omit(DEG_glia)
# nrDEG_glia=DEG_glia[DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.6,]
# choose_gene_neuron=rownames(nrDEG_glia)
# choose_matrix=group_D3A3[choose_gene_neuron,]
# DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
# up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0),]
# down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0),]
# setwd("~/Documents/Virginia Tech/Project_3_MIA mouse/differential analysis/results/DEGs/ABCD/union_FC_0.6")
# all_glia = group_D3A3[rownames(DEG_glia),]
# all_matrix = data.frame(DEG_glia,all_glia)
# write.csv(all_matrix,'All_D3A3.csv')
# write.csv(DEG_matrix,'Differential_D3A3.csv')
