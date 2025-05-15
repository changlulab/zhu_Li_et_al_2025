###########preparation##############
library(dplyr)
library(DESeq2)
library(umap)
library(ggplot2)
rm(list = ls())
options(stringsAsFactors = F)
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/")
factors = read.csv('cofactors_promoter.csv')
row.names(factors) = factors[,1]
factors = factors[,-1]
factors$Peaks = round(factors$Peaks,2)
factors$sequencing_depth = round(factors$sequencing_depth,2)
factors$unique_reads = round(factors$unique_reads,2)
##########################ABCD####################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/count_files/")
raw_def = read.csv('promoter_abcd_count.csv')
rownames(raw_def) = raw_def[,4]
def = raw_def[,5:ncol(raw_def)]
peak_info = raw_def[,1:4]
####strong correlation: cor>0.6################
###############################B vs. A_1###############################################################
meta = bind_rows(factors[29:40,1:8],factors[1:12,1:8])
meta$condition_num = c(rep(1,12),rep(2,12))
group_BA1 = data.frame(def[,29:40],def[,1:12])
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
cor_con_peaks = cor(meta$condition_num,meta$Peaks,method = 'pearson')
cor_con_NSC = cor(meta$condition_num,meta$NSC,method = 'pearson')
cor_con_RSC = cor(meta$condition_num,meta$RSC,method = 'pearson')
cor_con_FRiP = cor(meta$condition_num,meta$FRiP,method = 'pearson')
cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
cor_con_uni_reads = cor(meta$condition_num,meta$unique_reads,method = 'pearson')
#  
# NSC -UMAP
corr_uni_reads_1 = cor(umaps$X1,meta$unique_reads,method = 'pearson')
corr_uni_reads_2 = cor(umaps$X2,meta$unique_reads,method = 'pearson')
corr_uni_reads_3 = cor(umaps$X3,meta$unique_reads,method = 'pearson')
corr_uni_reads_4 = cor(umaps$X4,meta$unique_reads,method = 'pearson')
corr_uni_reads_5 = cor(umaps$X5,meta$unique_reads,method = 'pearson')
corr_uni_reads_6 = cor(umaps$X6,meta$unique_reads,method = 'pearson')

corr_peaks_1 = cor(umaps$X1,meta$Peaks,method = 'pearson')
corr_peaks_2 = cor(umaps$X2,meta$Peaks,method = 'pearson')
corr_peaks_3 = cor(umaps$X3,meta$Peaks,method = 'pearson')
corr_peaks_4 = cor(umaps$X4,meta$Peaks,method = 'pearson')
corr_peaks_5 = cor(umaps$X5,meta$Peaks,method = 'pearson')
corr_peaks_6 = cor(umaps$X6,meta$Peaks,method = 'pearson')

corr_NSC_1 = cor(umaps$X1,meta$NSC,method = 'pearson')
corr_NSC_2 = cor(umaps$X2,meta$NSC,method = 'pearson')
corr_NSC_3 = cor(umaps$X3,meta$NSC,method = 'pearson')
corr_NSC_4 = cor(umaps$X4,meta$NSC,method = 'pearson')
corr_NSC_5 = cor(umaps$X5,meta$NSC,method = 'pearson')
corr_NSC_6 = cor(umaps$X6,meta$NSC,method = 'pearson')

corr_RSC_1 = cor(umaps$X1,meta$RSC,method = 'pearson')
corr_RSC_2 = cor(umaps$X2,meta$RSC,method = 'pearson')
corr_RSC_3 = cor(umaps$X3,meta$RSC,method = 'pearson')
corr_RSC_4 = cor(umaps$X4,meta$RSC,method = 'pearson')
corr_RSC_5 = cor(umaps$X5,meta$RSC,method = 'pearson')
corr_RSC_6 = cor(umaps$X6,meta$RSC,method = 'pearson')

corr_FRiP_1 = cor(umaps$X1,meta$FRiP,method = 'pearson')
corr_FRiP_2 = cor(umaps$X2,meta$FRiP,method = 'pearson')
corr_FRiP_3 = cor(umaps$X3,meta$FRiP,method = 'pearson')
corr_FRiP_4 = cor(umaps$X4,meta$FRiP,method = 'pearson')
corr_FRiP_5 = cor(umaps$X5,meta$FRiP,method = 'pearson')
corr_FRiP_6 = cor(umaps$X6,meta$FRiP,method = 'pearson')

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
group_list2 = meta$unique_reads
group_list3 = meta$Peaks
group_list4 = meta$NSC
group_list5 = meta$RSC
group_list6 = meta$FRiP
group_list7 = meta$sequencing_depth
group_list8 = meta$batch

colData=data.frame(row.names = colnames(group_BA1),
                   condition=group_list1,
                   reads=group_list2,
                   peaks=group_list3,
                   NSC=group_list4,
                   RSC = group_list5,
                   FRiP = group_list6,
                   sequencing_depth = group_list7,
                   batch = group_list8)

colData$condition = factor(colData$condition,levels = c('B','A'))
colData$reads =as.numeric(colData$reads)
colData$peaks=as.numeric(colData$peaks)
colData$NSC=as.numeric(colData$NSC)
colData$RSC=as.numeric(colData$RSC)
colData$FRiP=as.numeric(colData$FRiP)
colData$batch=as.factor(colData$batch)
colData$sequencing_depth=as.numeric(colData$sequencing_depth)

dds=DESeqDataSetFromMatrix(countData = group_BA1,
                           colData = colData,
                           design = ~ condition + FRiP + NSC + RSC)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","B","A"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.6),]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_BA1[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0.5),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < -0.5),]
diff_peak = peak_info[rownames(DEG_matrix),]
new_DEG_matrix = bind_cols(diff_peak,DEG_matrix)
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/BA1")
all_glia = group_BA1[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_BA1.csv')
write.csv(new_DEG_matrix,'Differential_BA1.csv')
write.table(new_DEG_matrix,'Differential_BA1.bed')
###############################B vs. A_2###############################################################
meta = bind_rows(factors[29:40,1:8],factors[13:18,1:8])
meta$condition_num = c(rep(1,12),rep(2,6))
group_BA2 = data.frame(def[,29:40],def[,13:18])
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
cor_con_peaks = cor(meta$condition_num,meta$Peaks,method = 'pearson')
cor_con_NSC = cor(meta$condition_num,meta$NSC,method = 'pearson')
cor_con_RSC = cor(meta$condition_num,meta$RSC,method = 'pearson')
cor_con_FRiP = cor(meta$condition_num,meta$FRiP,method = 'pearson')
cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
cor_con_uni_reads = cor(meta$condition_num,meta$unique_reads,method = 'pearson')
# peaks RSC reads 
# batch peaks -UMAP
corr_uni_reads_1 = cor(umaps$X1,meta$unique_reads,method = 'pearson')
corr_uni_reads_2 = cor(umaps$X2,meta$unique_reads,method = 'pearson')
corr_uni_reads_3 = cor(umaps$X3,meta$unique_reads,method = 'pearson')
corr_uni_reads_4 = cor(umaps$X4,meta$unique_reads,method = 'pearson')
corr_uni_reads_5 = cor(umaps$X5,meta$unique_reads,method = 'pearson')
corr_uni_reads_6 = cor(umaps$X6,meta$unique_reads,method = 'pearson')

corr_peaks_1 = cor(umaps$X1,meta$Peaks,method = 'pearson')
corr_peaks_2 = cor(umaps$X2,meta$Peaks,method = 'pearson')
corr_peaks_3 = cor(umaps$X3,meta$Peaks,method = 'pearson')
corr_peaks_4 = cor(umaps$X4,meta$Peaks,method = 'pearson')
corr_peaks_5 = cor(umaps$X5,meta$Peaks,method = 'pearson')
corr_peaks_6 = cor(umaps$X6,meta$Peaks,method = 'pearson')

corr_NSC_1 = cor(umaps$X1,meta$NSC,method = 'pearson')
corr_NSC_2 = cor(umaps$X2,meta$NSC,method = 'pearson')
corr_NSC_3 = cor(umaps$X3,meta$NSC,method = 'pearson')
corr_NSC_4 = cor(umaps$X4,meta$NSC,method = 'pearson')
corr_NSC_5 = cor(umaps$X5,meta$NSC,method = 'pearson')
corr_NSC_6 = cor(umaps$X6,meta$NSC,method = 'pearson')

corr_RSC_1 = cor(umaps$X1,meta$RSC,method = 'pearson')
corr_RSC_2 = cor(umaps$X2,meta$RSC,method = 'pearson')
corr_RSC_3 = cor(umaps$X3,meta$RSC,method = 'pearson')
corr_RSC_4 = cor(umaps$X4,meta$RSC,method = 'pearson')
corr_RSC_5 = cor(umaps$X5,meta$RSC,method = 'pearson')
corr_RSC_6 = cor(umaps$X6,meta$RSC,method = 'pearson')

corr_FRiP_1 = cor(umaps$X1,meta$FRiP,method = 'pearson')
corr_FRiP_2 = cor(umaps$X2,meta$FRiP,method = 'pearson')
corr_FRiP_3 = cor(umaps$X3,meta$FRiP,method = 'pearson')
corr_FRiP_4 = cor(umaps$X4,meta$FRiP,method = 'pearson')
corr_FRiP_5 = cor(umaps$X5,meta$FRiP,method = 'pearson')
corr_FRiP_6 = cor(umaps$X6,meta$FRiP,method = 'pearson')

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
group_list2 = meta$unique_reads
group_list3 = meta$Peaks
group_list4 = meta$NSC
group_list5 = meta$RSC
group_list6 = meta$FRiP
group_list7 = meta$sequencing_depth
group_list8 = meta$batch

colData=data.frame(row.names = colnames(group_BA2),
                   condition=group_list1,
                   reads=group_list2,
                   peaks=group_list3,
                   NSC=group_list4,
                   RSC = group_list5,
                   FRiP = group_list6,
                   sequencing_depth = group_list7,
                   batch = group_list8)

colData$condition = factor(colData$condition,levels = c('B','A'))
colData$reads =as.numeric(colData$reads)
colData$peaks=as.numeric(colData$peaks)
colData$NSC=as.numeric(colData$NSC)
colData$RSC=as.numeric(colData$RSC)
colData$FRiP=as.numeric(colData$FRiP)
colData$sequencing_depth=as.numeric(colData$sequencing_depth)
colData$batch=as.factor(colData$batch)

dds=DESeqDataSetFromMatrix(countData = group_BA2,
                           colData = colData,
                           design = ~ condition + peaks + reads)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","B","A"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.6),]
# nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.5),]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_BA2[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0 & nrDEG_glia$log2FoldChange > 0.5),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0 & nrDEG_glia$log2FoldChange < -0.5),]
diff_peak = peak_info[rownames(DEG_matrix),]
new_DEG_matrix = bind_cols(diff_peak,DEG_matrix)
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/BA2")
all_glia = group_BA2[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_BA2.csv')
write.csv(new_DEG_matrix,'Differential_BA2.csv')
###############################B vs. A_3###############################################################
meta = bind_rows(factors[29:40,1:8],factors[19:28,1:8])
meta$condition_num = c(rep(1,12),rep(2,10))
group_BA3 = data.frame(def[,29:40],def[,19:28])
miss_1 <- c()
for (i in 1:nrow(group_BA3)) {
  if(length(which(group_BA3[i,1:12]<10)) > 0.5*ncol(group_BA3[,1:12])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(group_BA3)) {
  if(length(which(group_BA3[i,13:22]<10)) > 0.5*ncol(group_BA3[,13:22])) miss_2 <- append(miss_2,i)
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
cor_con_peaks = cor(meta$condition_num,meta$Peaks,method = 'pearson')
cor_con_NSC = cor(meta$condition_num,meta$NSC,method = 'pearson')
cor_con_RSC = cor(meta$condition_num,meta$RSC,method = 'pearson')
cor_con_FRiP = cor(meta$condition_num,meta$FRiP,method = 'pearson')
cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
cor_con_uni_reads = cor(meta$condition_num,meta$unique_reads,method = 'pearson')
# NSC peaks reads 
# -UMAP
corr_uni_reads_1 = cor(umaps$X1,meta$unique_reads,method = 'pearson')
corr_uni_reads_2 = cor(umaps$X2,meta$unique_reads,method = 'pearson')
corr_uni_reads_3 = cor(umaps$X3,meta$unique_reads,method = 'pearson')
corr_uni_reads_4 = cor(umaps$X4,meta$unique_reads,method = 'pearson')
corr_uni_reads_5 = cor(umaps$X5,meta$unique_reads,method = 'pearson')
corr_uni_reads_6 = cor(umaps$X6,meta$unique_reads,method = 'pearson')

corr_peaks_1 = cor(umaps$X1,meta$Peaks,method = 'pearson')
corr_peaks_2 = cor(umaps$X2,meta$Peaks,method = 'pearson')
corr_peaks_3 = cor(umaps$X3,meta$Peaks,method = 'pearson')
corr_peaks_4 = cor(umaps$X4,meta$Peaks,method = 'pearson')
corr_peaks_5 = cor(umaps$X5,meta$Peaks,method = 'pearson')
corr_peaks_6 = cor(umaps$X6,meta$Peaks,method = 'pearson')

corr_NSC_1 = cor(umaps$X1,meta$NSC,method = 'pearson')
corr_NSC_2 = cor(umaps$X2,meta$NSC,method = 'pearson')
corr_NSC_3 = cor(umaps$X3,meta$NSC,method = 'pearson')
corr_NSC_4 = cor(umaps$X4,meta$NSC,method = 'pearson')
corr_NSC_5 = cor(umaps$X5,meta$NSC,method = 'pearson')
corr_NSC_6 = cor(umaps$X6,meta$NSC,method = 'pearson')

corr_RSC_1 = cor(umaps$X1,meta$RSC,method = 'pearson')
corr_RSC_2 = cor(umaps$X2,meta$RSC,method = 'pearson')
corr_RSC_3 = cor(umaps$X3,meta$RSC,method = 'pearson')
corr_RSC_4 = cor(umaps$X4,meta$RSC,method = 'pearson')
corr_RSC_5 = cor(umaps$X5,meta$RSC,method = 'pearson')
corr_RSC_6 = cor(umaps$X6,meta$RSC,method = 'pearson')

corr_FRiP_1 = cor(umaps$X1,meta$FRiP,method = 'pearson')
corr_FRiP_2 = cor(umaps$X2,meta$FRiP,method = 'pearson')
corr_FRiP_3 = cor(umaps$X3,meta$FRiP,method = 'pearson')
corr_FRiP_4 = cor(umaps$X4,meta$FRiP,method = 'pearson')
corr_FRiP_5 = cor(umaps$X5,meta$FRiP,method = 'pearson')
corr_FRiP_6 = cor(umaps$X6,meta$FRiP,method = 'pearson')

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
group_list2 = meta$unique_reads
group_list3 = meta$Peaks
group_list4 = meta$NSC
group_list5 = meta$RSC
group_list6 = meta$FRiP
group_list7 = meta$sequencing_depth
group_list8 = meta$batch

colData=data.frame(row.names = colnames(group_BA3),
                   condition=group_list1,
                   reads=group_list2,
                   peaks=group_list3,
                   NSC=group_list4,
                   RSC = group_list5,
                   FRiP = group_list6,
                   sequencing_depth = group_list7,
                   batch = group_list8)

colData$condition = factor(colData$condition,levels = c('B','A'))
colData$reads =as.numeric(colData$reads)
colData$peaks=as.numeric(colData$peaks)
colData$NSC=as.numeric(colData$NSC)
colData$RSC=as.numeric(colData$RSC)
colData$batch=as.factor(colData$batch)
colData$FRiP=as.numeric(colData$FRiP)
colData$sequencing_depth=as.numeric(colData$sequencing_depth)

dds=DESeqDataSetFromMatrix(countData = group_BA3,
                           colData = colData,
                           design = ~ condition + NSC + sequencing_depth + reads)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","B","A"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.6),]
# nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.5),]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_BA3[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0 & nrDEG_glia$log2FoldChange > 0.5),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0 & nrDEG_glia$log2FoldChange < -0.5),]
diff_peak = peak_info[rownames(DEG_matrix),]
new_DEG_matrix = bind_cols(diff_peak,DEG_matrix)
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/BA3")
all_glia = group_BA3[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_BA3.csv')
write.csv(new_DEG_matrix,'Differential_BA3.csv')
##############################D_1 vs. C###############################################################
meta = bind_rows(factors[53:64,1:8],factors[41:52,1:8])
meta$condition_num = c(rep(1,12),rep(2,12))
group_D1C = data.frame(def[,53:64],def[,41:52])
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
cor_con_peaks = cor(meta$condition_num,meta$Peaks,method = 'pearson')
cor_con_NSC = cor(meta$condition_num,meta$NSC,method = 'pearson')
cor_con_RSC = cor(meta$condition_num,meta$RSC,method = 'pearson')
cor_con_FRiP = cor(meta$condition_num,meta$FRiP,method = 'pearson')
cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
cor_con_uni_reads = cor(meta$condition_num,meta$unique_reads,method = 'pearson')
#
# FRiP NSC -UMAP
corr_uni_reads_1 = cor(umaps$X1,meta$unique_reads,method = 'pearson')
corr_uni_reads_2 = cor(umaps$X2,meta$unique_reads,method = 'pearson')
corr_uni_reads_3 = cor(umaps$X3,meta$unique_reads,method = 'pearson')
corr_uni_reads_4 = cor(umaps$X4,meta$unique_reads,method = 'pearson')
corr_uni_reads_5 = cor(umaps$X5,meta$unique_reads,method = 'pearson')
corr_uni_reads_6 = cor(umaps$X6,meta$unique_reads,method = 'pearson')

corr_peaks_1 = cor(umaps$X1,meta$Peaks,method = 'pearson')
corr_peaks_2 = cor(umaps$X2,meta$Peaks,method = 'pearson')
corr_peaks_3 = cor(umaps$X3,meta$Peaks,method = 'pearson')
corr_peaks_4 = cor(umaps$X4,meta$Peaks,method = 'pearson')
corr_peaks_5 = cor(umaps$X5,meta$Peaks,method = 'pearson')
corr_peaks_6 = cor(umaps$X6,meta$Peaks,method = 'pearson')

corr_NSC_1 = cor(umaps$X1,meta$NSC,method = 'pearson')
corr_NSC_2 = cor(umaps$X2,meta$NSC,method = 'pearson')
corr_NSC_3 = cor(umaps$X3,meta$NSC,method = 'pearson')
corr_NSC_4 = cor(umaps$X4,meta$NSC,method = 'pearson')
corr_NSC_5 = cor(umaps$X5,meta$NSC,method = 'pearson')
corr_NSC_6 = cor(umaps$X6,meta$NSC,method = 'pearson')

corr_RSC_1 = cor(umaps$X1,meta$RSC,method = 'pearson')
corr_RSC_2 = cor(umaps$X2,meta$RSC,method = 'pearson')
corr_RSC_3 = cor(umaps$X3,meta$RSC,method = 'pearson')
corr_RSC_4 = cor(umaps$X4,meta$RSC,method = 'pearson')
corr_RSC_5 = cor(umaps$X5,meta$RSC,method = 'pearson')
corr_RSC_6 = cor(umaps$X6,meta$RSC,method = 'pearson')

corr_FRiP_1 = cor(umaps$X1,meta$FRiP,method = 'pearson')
corr_FRiP_2 = cor(umaps$X2,meta$FRiP,method = 'pearson')
corr_FRiP_3 = cor(umaps$X3,meta$FRiP,method = 'pearson')
corr_FRiP_4 = cor(umaps$X4,meta$FRiP,method = 'pearson')
corr_FRiP_5 = cor(umaps$X5,meta$FRiP,method = 'pearson')
corr_FRiP_6 = cor(umaps$X6,meta$FRiP,method = 'pearson')

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
group_list2 = meta$unique_reads
group_list3 = meta$Peaks
group_list4 = meta$NSC
group_list5 = meta$RSC
group_list6 = meta$FRiP
group_list7 = meta$sequencing_depth
group_list8 = meta$batch

colData=data.frame(row.names = colnames(group_D1C),
                   condition=group_list1,
                   reads=group_list2,
                   peaks=group_list3,
                   NSC=group_list4,
                   RSC = group_list5,
                   FRiP = group_list6,
                   sequencing_depth = group_list7,
                   batch = group_list8)

colData$condition = factor(colData$condition,levels = c('D','C'))
colData$reads =as.numeric(colData$reads)
colData$peaks=as.numeric(colData$peaks)
colData$NSC=as.numeric(colData$NSC)
colData$RSC=as.numeric(colData$RSC)
colData$FRiP=as.numeric(colData$FRiP)
colData$sequencing_depth=as.numeric(colData$sequencing_depth)
colData$batch=as.numeric(colData$batch)

dds=DESeqDataSetFromMatrix(countData = group_D1C,
                           colData = colData,
                           design = ~ condition + RSC + NSC + FRiP)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","D","C"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.6),]
# nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.5),]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_D1C[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0 & nrDEG_glia$log2FoldChange > 0.5),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0 & nrDEG_glia$log2FoldChange < -0.5),]
diff_peak = peak_info[rownames(DEG_matrix),]
new_DEG_matrix = bind_cols(diff_peak,DEG_matrix)
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/D1C")
all_glia = group_D1C[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_D1C.csv')
write.csv(new_DEG_matrix,'Differential_D1C.csv')
###############################D_2 vs. C###############################################################
meta = bind_rows(factors[65:70,1:8],factors[41:52,1:8])
meta$condition_num = c(rep(1,6),rep(2,12))
group_D2C = data.frame(def[,65:70],def[,41:52])
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
cor_con_peaks = cor(meta$condition_num,meta$Peaks,method = 'pearson')
cor_con_NSC = cor(meta$condition_num,meta$NSC,method = 'pearson')
cor_con_RSC = cor(meta$condition_num,meta$RSC,method = 'pearson')
cor_con_FRiP = cor(meta$condition_num,meta$FRiP,method = 'pearson')
cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
cor_con_uni_reads = cor(meta$condition_num,meta$unique_reads,method = 'pearson')
#
# peaks -UMAP
corr_uni_reads_1 = cor(umaps$X1,meta$unique_reads,method = 'pearson')
corr_uni_reads_2 = cor(umaps$X2,meta$unique_reads,method = 'pearson')
corr_uni_reads_3 = cor(umaps$X3,meta$unique_reads,method = 'pearson')
corr_uni_reads_4 = cor(umaps$X4,meta$unique_reads,method = 'pearson')
corr_uni_reads_5 = cor(umaps$X5,meta$unique_reads,method = 'pearson')
corr_uni_reads_6 = cor(umaps$X6,meta$unique_reads,method = 'pearson')

corr_peaks_1 = cor(umaps$X1,meta$Peaks,method = 'pearson')
corr_peaks_2 = cor(umaps$X2,meta$Peaks,method = 'pearson')
corr_peaks_3 = cor(umaps$X3,meta$Peaks,method = 'pearson')
corr_peaks_4 = cor(umaps$X4,meta$Peaks,method = 'pearson')
corr_peaks_5 = cor(umaps$X5,meta$Peaks,method = 'pearson')
corr_peaks_6 = cor(umaps$X6,meta$Peaks,method = 'pearson')

corr_NSC_1 = cor(umaps$X1,meta$NSC,method = 'pearson')
corr_NSC_2 = cor(umaps$X2,meta$NSC,method = 'pearson')
corr_NSC_3 = cor(umaps$X3,meta$NSC,method = 'pearson')
corr_NSC_4 = cor(umaps$X4,meta$NSC,method = 'pearson')
corr_NSC_5 = cor(umaps$X5,meta$NSC,method = 'pearson')
corr_NSC_6 = cor(umaps$X6,meta$NSC,method = 'pearson')

corr_RSC_1 = cor(umaps$X1,meta$RSC,method = 'pearson')
corr_RSC_2 = cor(umaps$X2,meta$RSC,method = 'pearson')
corr_RSC_3 = cor(umaps$X3,meta$RSC,method = 'pearson')
corr_RSC_4 = cor(umaps$X4,meta$RSC,method = 'pearson')
corr_RSC_5 = cor(umaps$X5,meta$RSC,method = 'pearson')
corr_RSC_6 = cor(umaps$X6,meta$RSC,method = 'pearson')

corr_FRiP_1 = cor(umaps$X1,meta$FRiP,method = 'pearson')
corr_FRiP_2 = cor(umaps$X2,meta$FRiP,method = 'pearson')
corr_FRiP_3 = cor(umaps$X3,meta$FRiP,method = 'pearson')
corr_FRiP_4 = cor(umaps$X4,meta$FRiP,method = 'pearson')
corr_FRiP_5 = cor(umaps$X5,meta$FRiP,method = 'pearson')
corr_FRiP_6 = cor(umaps$X6,meta$FRiP,method = 'pearson')

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
group_list2 = meta$unique_reads
group_list3 = meta$Peaks
group_list4 = meta$NSC
group_list5 = meta$RSC
group_list6 = meta$FRiP
group_list7 = meta$sequencing_depth
group_list8 = meta$batch

colData=data.frame(row.names = colnames(group_D2C),
                   condition=group_list1,
                   reads=group_list2,
                   peaks=group_list3,
                   NSC=group_list4,
                   RSC = group_list5,
                   FRiP = group_list6,
                   sequencing_depth = group_list7,
                   batch = group_list8)

colData$condition = factor(colData$condition,levels = c('D','C'))
colData$reads =as.numeric(colData$reads)
colData$peaks=as.numeric(colData$peaks)
colData$NSC=as.numeric(colData$NSC)
colData$RSC=as.numeric(colData$RSC)
colData$FRiP=as.numeric(colData$FRiP)
colData$sequencing_depth=as.numeric(colData$sequencing_depth)
colData$batch=as.factor(colData$batch)

dds=DESeqDataSetFromMatrix(countData = group_D2C,
                           colData = colData,
                           design = ~ condition + peaks + NSC)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","D","C"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.6),]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_D2C[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0 & nrDEG_glia$log2FoldChange > 0.5),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0 & nrDEG_glia$log2FoldChange < -0.5),]
diff_peak = peak_info[rownames(DEG_matrix),]
new_DEG_matrix = bind_cols(diff_peak,DEG_matrix)
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/D2C")
all_glia = group_D2C[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_D2C.csv')
write.csv(new_DEG_matrix,'Differential_D2C.csv')
###############################D_3 vs. C###############################################################
meta = bind_rows(factors[71:82,1:8],factors[41:52,1:8])
meta$condition_num = c(rep(1,12),rep(2,12))
group_D3C = data.frame(def[,71:82],def[,41:52])
miss_1 <- c()
for (i in 1:nrow(group_D3C)) {
  if(length(which(group_D3C[i,1:12]<10)) > 0.5*ncol(group_D3C[,1:12])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(group_D3C)) {
  if(length(which(group_D3C[i,13:24]<10)) > 0.5*ncol(group_D3C[,13:24])) miss_2 <- append(miss_2,i)
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
cor_con_peaks = cor(meta$condition_num,meta$Peaks,method = 'pearson')
cor_con_NSC = cor(meta$condition_num,meta$NSC,method = 'pearson')
cor_con_RSC = cor(meta$condition_num,meta$RSC,method = 'pearson')
cor_con_FRiP = cor(meta$condition_num,meta$FRiP,method = 'pearson')
cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
cor_con_uni_reads = cor(meta$condition_num,meta$unique_reads,method = 'pearson')
# NSC peaks
# NSC peaks -UMAP
corr_uni_reads_1 = cor(umaps$X1,meta$unique_reads,method = 'pearson')
corr_uni_reads_2 = cor(umaps$X2,meta$unique_reads,method = 'pearson')
corr_uni_reads_3 = cor(umaps$X3,meta$unique_reads,method = 'pearson')
corr_uni_reads_4 = cor(umaps$X4,meta$unique_reads,method = 'pearson')
corr_uni_reads_5 = cor(umaps$X5,meta$unique_reads,method = 'pearson')
corr_uni_reads_6 = cor(umaps$X6,meta$unique_reads,method = 'pearson')

corr_peaks_1 = cor(umaps$X1,meta$Peaks,method = 'pearson')
corr_peaks_2 = cor(umaps$X2,meta$Peaks,method = 'pearson')
corr_peaks_3 = cor(umaps$X3,meta$Peaks,method = 'pearson')
corr_peaks_4 = cor(umaps$X4,meta$Peaks,method = 'pearson')
corr_peaks_5 = cor(umaps$X5,meta$Peaks,method = 'pearson')
corr_peaks_6 = cor(umaps$X6,meta$Peaks,method = 'pearson')

corr_NSC_1 = cor(umaps$X1,meta$NSC,method = 'pearson')
corr_NSC_2 = cor(umaps$X2,meta$NSC,method = 'pearson')
corr_NSC_3 = cor(umaps$X3,meta$NSC,method = 'pearson')
corr_NSC_4 = cor(umaps$X4,meta$NSC,method = 'pearson')
corr_NSC_5 = cor(umaps$X5,meta$NSC,method = 'pearson')
corr_NSC_6 = cor(umaps$X6,meta$NSC,method = 'pearson')

corr_RSC_1 = cor(umaps$X1,meta$RSC,method = 'pearson')
corr_RSC_2 = cor(umaps$X2,meta$RSC,method = 'pearson')
corr_RSC_3 = cor(umaps$X3,meta$RSC,method = 'pearson')
corr_RSC_4 = cor(umaps$X4,meta$RSC,method = 'pearson')
corr_RSC_5 = cor(umaps$X5,meta$RSC,method = 'pearson')
corr_RSC_6 = cor(umaps$X6,meta$RSC,method = 'pearson')

corr_FRiP_1 = cor(umaps$X1,meta$FRiP,method = 'pearson')
corr_FRiP_2 = cor(umaps$X2,meta$FRiP,method = 'pearson')
corr_FRiP_3 = cor(umaps$X3,meta$FRiP,method = 'pearson')
corr_FRiP_4 = cor(umaps$X4,meta$FRiP,method = 'pearson')
corr_FRiP_5 = cor(umaps$X5,meta$FRiP,method = 'pearson')
corr_FRiP_6 = cor(umaps$X6,meta$FRiP,method = 'pearson')

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
group_list2 = meta$unique_reads
group_list3 = meta$Peaks
group_list4 = meta$NSC
group_list5 = meta$RSC
group_list6 = meta$FRiP
group_list7 = meta$sequencing_depth
group_list8 = meta$batch

colData=data.frame(row.names = colnames(group_D3C),
                   condition=group_list1,
                   reads=group_list2,
                   peaks=group_list3,
                   NSC=group_list4,
                   RSC = group_list5,
                   FRiP = group_list6,
                   sequencing_depth = group_list7,
                   batch = group_list8)

colData$condition = factor(colData$condition,levels = c('D','C'))
colData$reads =as.numeric(colData$reads)
colData$peaks=as.numeric(colData$peaks)
colData$NSC=as.numeric(colData$NSC)
colData$RSC=as.numeric(colData$RSC)
colData$FRiP=as.numeric(colData$FRiP)
colData$sequencing_depth=as.numeric(colData$sequencing_depth)
colData$batch=as.numeric(colData$batch)

dds=DESeqDataSetFromMatrix(countData = group_D3C,
                           colData = colData,
                           design = ~ condition + NSC + batch)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","D","C"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) >0.6),]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_D3C[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0 & nrDEG_glia$log2FoldChange > 0.5),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0 & nrDEG_glia$log2FoldChange < -0.5),]
diff_peak = peak_info[rownames(DEG_matrix),]
new_DEG_matrix = bind_cols(diff_peak,DEG_matrix)
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/D3C")
all_glia = group_D3C[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_D3C.csv')
write.csv(new_DEG_matrix,'Differential_D3C.csv')
# ###############################C vs. A_1###############################################################
meta = bind_rows(factors[41:52,1:8],factors[1:12,1:8])
meta$condition_num = c(rep(1,12),rep(2,12))
group_CA1 = data.frame(def[,41:52],def[,1:12])
miss_1 <- c()
for (i in 1:nrow(group_CA1)) {
  if(length(which(group_CA1[i,1:12]<12)) > 0.5*ncol(group_CA1[,1:12])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(group_CA1)) {
  if(length(which(group_CA1[i,13:24]<10)) > 0.5*ncol(group_CA1[,13:24])) miss_2 <- append(miss_2,i)
}
miss = union(miss_1,miss_2)
group_CA1 = group_CA1[-miss,]
t_group_CA1 = t(group_CA1)
umap_results <- umap(t_group_CA1,n_neighbors=11,min_dist=0.1,n_components = 6)
umaps = data.frame(umap_results[["layout"]])
umaps$condition = as.factor(meta$condition)
ggplot(umaps,mapping = aes(x=X1,y=X2,col=condition)) +
  geom_point(aes(colour = condition))

cor_con_batch = cor(meta$condition_num,meta$batch,method = 'pearson')
cor_con_peaks = cor(meta$condition_num,meta$Peaks,method = 'pearson')
cor_con_NSC = cor(meta$condition_num,meta$NSC,method = 'pearson')
cor_con_RSC = cor(meta$condition_num,meta$RSC,method = 'pearson')
cor_con_FRiP = cor(meta$condition_num,meta$FRiP,method = 'pearson')
cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
cor_con_uni_reads = cor(meta$condition_num,meta$unique_reads,method = 'pearson')
#
# FRiP NSC  -UMAP
corr_uni_reads_1 = cor(umaps$X1,meta$unique_reads,method = 'pearson')
corr_uni_reads_2 = cor(umaps$X2,meta$unique_reads,method = 'pearson')
corr_uni_reads_3 = cor(umaps$X3,meta$unique_reads,method = 'pearson')
corr_uni_reads_4 = cor(umaps$X4,meta$unique_reads,method = 'pearson')
corr_uni_reads_5 = cor(umaps$X5,meta$unique_reads,method = 'pearson')
corr_uni_reads_6 = cor(umaps$X6,meta$unique_reads,method = 'pearson')

corr_peaks_1 = cor(umaps$X1,meta$Peaks,method = 'pearson')
corr_peaks_2 = cor(umaps$X2,meta$Peaks,method = 'pearson')
corr_peaks_3 = cor(umaps$X3,meta$Peaks,method = 'pearson')
corr_peaks_4 = cor(umaps$X4,meta$Peaks,method = 'pearson')
corr_peaks_5 = cor(umaps$X5,meta$Peaks,method = 'pearson')
corr_peaks_6 = cor(umaps$X6,meta$Peaks,method = 'pearson')

corr_NSC_1 = cor(umaps$X1,meta$NSC,method = 'pearson')
corr_NSC_2 = cor(umaps$X2,meta$NSC,method = 'pearson')
corr_NSC_3 = cor(umaps$X3,meta$NSC,method = 'pearson')
corr_NSC_4 = cor(umaps$X4,meta$NSC,method = 'pearson')
corr_NSC_5 = cor(umaps$X5,meta$NSC,method = 'pearson')
corr_NSC_6 = cor(umaps$X6,meta$NSC,method = 'pearson')

corr_RSC_1 = cor(umaps$X1,meta$RSC,method = 'pearson')
corr_RSC_2 = cor(umaps$X2,meta$RSC,method = 'pearson')
corr_RSC_3 = cor(umaps$X3,meta$RSC,method = 'pearson')
corr_RSC_4 = cor(umaps$X4,meta$RSC,method = 'pearson')
corr_RSC_5 = cor(umaps$X5,meta$RSC,method = 'pearson')
corr_RSC_6 = cor(umaps$X6,meta$RSC,method = 'pearson')

corr_FRiP_1 = cor(umaps$X1,meta$FRiP,method = 'pearson')
corr_FRiP_2 = cor(umaps$X2,meta$FRiP,method = 'pearson')
corr_FRiP_3 = cor(umaps$X3,meta$FRiP,method = 'pearson')
corr_FRiP_4 = cor(umaps$X4,meta$FRiP,method = 'pearson')
corr_FRiP_5 = cor(umaps$X5,meta$FRiP,method = 'pearson')
corr_FRiP_6 = cor(umaps$X6,meta$FRiP,method = 'pearson')

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
group_list2 = meta$unique_reads
group_list3 = meta$Peaks
group_list4 = meta$NSC
group_list5 = meta$RSC
group_list6 = meta$FRiP
group_list7 = meta$sequencing_depth
group_list8 = meta$batch

colData=data.frame(row.names = colnames(group_CA1),
                   condition=group_list1,
                   reads=group_list2,
                   peaks=group_list3,
                   NSC=group_list4,
                   RSC = group_list5,
                   FRiP = group_list6,
                   sequencing_depth = group_list7,
                   batch = group_list8)

colData$condition = factor(colData$condition,levels = c('C','A'))
colData$reads =as.numeric(colData$reads)
colData$peaks=as.numeric(colData$peaks)
colData$NSC=as.numeric(colData$NSC)
colData$RSC=as.numeric(colData$RSC)
colData$FRiP=as.numeric(colData$FRiP)
colData$sequencing_depth=as.numeric(colData$sequencing_depth)
colData$batch=as.factor(colData$batch)

dds=DESeqDataSetFromMatrix(countData = group_CA1,
                           colData = colData,
                           design = ~ condition  + NSC + FRiP + RSC)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","C","A"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.6),]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_CA1[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0 & nrDEG_glia$log2FoldChange > 0.5),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0 & nrDEG_glia$log2FoldChange < -0.5),]
diff_peak = peak_info[rownames(DEG_matrix),]
new_DEG_matrix = bind_cols(diff_peak,DEG_matrix)
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/CA1")
all_glia = group_CA1[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_CA1.csv')
write.csv(new_DEG_matrix,'Differential_CA1.csv')
# ###############################C vs. A_2###############################################################
meta = bind_rows(factors[41:52,1:8],factors[13:18,1:8])
meta$condition_num = c(rep(1,12),rep(2,6))
group_CA2 = data.frame(def[,41:52],def[,13:18])
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
cor_con_peaks = cor(meta$condition_num,meta$Peaks,method = 'pearson')
cor_con_NSC = cor(meta$condition_num,meta$NSC,method = 'pearson')
cor_con_RSC = cor(meta$condition_num,meta$RSC,method = 'pearson')
cor_con_FRiP = cor(meta$condition_num,meta$FRiP,method = 'pearson')
cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
cor_con_uni_reads = cor(meta$condition_num,meta$unique_reads,method = 'pearson')
# peaks
# NSC -UMAP
corr_uni_reads_1 = cor(umaps$X1,meta$unique_reads,method = 'pearson')
corr_uni_reads_2 = cor(umaps$X2,meta$unique_reads,method = 'pearson')
corr_uni_reads_3 = cor(umaps$X3,meta$unique_reads,method = 'pearson')
corr_uni_reads_4 = cor(umaps$X4,meta$unique_reads,method = 'pearson')
corr_uni_reads_5 = cor(umaps$X5,meta$unique_reads,method = 'pearson')
corr_uni_reads_6 = cor(umaps$X6,meta$unique_reads,method = 'pearson')

corr_peaks_1 = cor(umaps$X1,meta$Peaks,method = 'pearson')
corr_peaks_2 = cor(umaps$X2,meta$Peaks,method = 'pearson')
corr_peaks_3 = cor(umaps$X3,meta$Peaks,method = 'pearson')
corr_peaks_4 = cor(umaps$X4,meta$Peaks,method = 'pearson')
corr_peaks_5 = cor(umaps$X5,meta$Peaks,method = 'pearson')
corr_peaks_6 = cor(umaps$X6,meta$Peaks,method = 'pearson')

corr_NSC_1 = cor(umaps$X1,meta$NSC,method = 'pearson')
corr_NSC_2 = cor(umaps$X2,meta$NSC,method = 'pearson')
corr_NSC_3 = cor(umaps$X3,meta$NSC,method = 'pearson')
corr_NSC_4 = cor(umaps$X4,meta$NSC,method = 'pearson')
corr_NSC_5 = cor(umaps$X5,meta$NSC,method = 'pearson')
corr_NSC_6 = cor(umaps$X6,meta$NSC,method = 'pearson')

corr_RSC_1 = cor(umaps$X1,meta$RSC,method = 'pearson')
corr_RSC_2 = cor(umaps$X2,meta$RSC,method = 'pearson')
corr_RSC_3 = cor(umaps$X3,meta$RSC,method = 'pearson')
corr_RSC_4 = cor(umaps$X4,meta$RSC,method = 'pearson')
corr_RSC_5 = cor(umaps$X5,meta$RSC,method = 'pearson')
corr_RSC_6 = cor(umaps$X6,meta$RSC,method = 'pearson')

corr_FRiP_1 = cor(umaps$X1,meta$FRiP,method = 'pearson')
corr_FRiP_2 = cor(umaps$X2,meta$FRiP,method = 'pearson')
corr_FRiP_3 = cor(umaps$X3,meta$FRiP,method = 'pearson')
corr_FRiP_4 = cor(umaps$X4,meta$FRiP,method = 'pearson')
corr_FRiP_5 = cor(umaps$X5,meta$FRiP,method = 'pearson')
corr_FRiP_6 = cor(umaps$X6,meta$FRiP,method = 'pearson')

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
group_list2 = meta$unique_reads
group_list3 = meta$Peaks
group_list4 = meta$NSC
group_list5 = meta$RSC
group_list6 = meta$FRiP
group_list7 = meta$sequencing_depth
group_list8 = meta$batch

colData=data.frame(row.names = colnames(group_CA2),
                   condition=group_list1,
                   reads=group_list2,
                   peaks=group_list3,
                   NSC=group_list4,
                   RSC = group_list5,
                   FRiP = group_list6,
                   sequencing_depth = group_list7,
                   batch = group_list8)

colData$condition = factor(colData$condition,levels = c('C','A'))
colData$reads =as.numeric(colData$reads)
colData$peaks=as.numeric(colData$peaks)
colData$NSC=as.numeric(colData$NSC)
colData$RSC=as.numeric(colData$RSC)
colData$FRiP=as.numeric(colData$FRiP)
colData$sequencing_depth=as.numeric(colData$sequencing_depth)
colData$batch=as.factor(colData$batch)

dds=DESeqDataSetFromMatrix(countData = group_CA2,
                           colData = colData,
                           design = ~ condition)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","C","A"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.6),]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_CA2[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0 & nrDEG_glia$log2FoldChange > 0.5),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0 & nrDEG_glia$log2FoldChange < -0.5),]
diff_peak = peak_info[rownames(DEG_matrix),]
new_DEG_matrix = bind_cols(diff_peak,DEG_matrix)
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/CA2")
all_glia = group_CA2[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_CA2.csv')
write.csv(new_DEG_matrix,'Differential_CA2.csv')
# ###############################C vs. A_3###############################################################
meta = bind_rows(factors[41:52,1:8],factors[19:28,1:8])
meta$condition_num = c(rep(1,12),rep(2,10))
group_CA3 = data.frame(def[,41:52],def[,19:28])
miss_1 <- c()
for (i in 1:nrow(group_CA3)) {
  if(length(which(group_CA3[i,1:12]<10)) > 0.5*ncol(group_CA3[,1:12])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(group_CA3)) {
  if(length(which(group_CA3[i,13:22]<10)) > 0.5*ncol(group_CA3[,12:22])) miss_2 <- append(miss_2,i)
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
cor_con_peaks = cor(meta$condition_num,meta$Peaks,method = 'pearson')
cor_con_NSC = cor(meta$condition_num,meta$NSC,method = 'pearson')
cor_con_RSC = cor(meta$condition_num,meta$RSC,method = 'pearson')
cor_con_FRiP = cor(meta$condition_num,meta$FRiP,method = 'pearson')
cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
cor_con_uni_reads = cor(meta$condition_num,meta$unique_reads,method = 'pearson')
# NSC peaks
#  -UMAP
corr_uni_reads_1 = cor(umaps$X1,meta$unique_reads,method = 'pearson')
corr_uni_reads_2 = cor(umaps$X2,meta$unique_reads,method = 'pearson')
corr_uni_reads_3 = cor(umaps$X3,meta$unique_reads,method = 'pearson')
corr_uni_reads_4 = cor(umaps$X4,meta$unique_reads,method = 'pearson')
corr_uni_reads_5 = cor(umaps$X5,meta$unique_reads,method = 'pearson')
corr_uni_reads_6 = cor(umaps$X6,meta$unique_reads,method = 'pearson')

corr_peaks_1 = cor(umaps$X1,meta$Peaks,method = 'pearson')
corr_peaks_2 = cor(umaps$X2,meta$Peaks,method = 'pearson')
corr_peaks_3 = cor(umaps$X3,meta$Peaks,method = 'pearson')
corr_peaks_4 = cor(umaps$X4,meta$Peaks,method = 'pearson')
corr_peaks_5 = cor(umaps$X5,meta$Peaks,method = 'pearson')
corr_peaks_6 = cor(umaps$X6,meta$Peaks,method = 'pearson')

corr_NSC_1 = cor(umaps$X1,meta$NSC,method = 'pearson')
corr_NSC_2 = cor(umaps$X2,meta$NSC,method = 'pearson')
corr_NSC_3 = cor(umaps$X3,meta$NSC,method = 'pearson')
corr_NSC_4 = cor(umaps$X4,meta$NSC,method = 'pearson')
corr_NSC_5 = cor(umaps$X5,meta$NSC,method = 'pearson')
corr_NSC_6 = cor(umaps$X6,meta$NSC,method = 'pearson')

corr_RSC_1 = cor(umaps$X1,meta$RSC,method = 'pearson')
corr_RSC_2 = cor(umaps$X2,meta$RSC,method = 'pearson')
corr_RSC_3 = cor(umaps$X3,meta$RSC,method = 'pearson')
corr_RSC_4 = cor(umaps$X4,meta$RSC,method = 'pearson')
corr_RSC_5 = cor(umaps$X5,meta$RSC,method = 'pearson')
corr_RSC_6 = cor(umaps$X6,meta$RSC,method = 'pearson')

corr_FRiP_1 = cor(umaps$X1,meta$FRiP,method = 'pearson')
corr_FRiP_2 = cor(umaps$X2,meta$FRiP,method = 'pearson')
corr_FRiP_3 = cor(umaps$X3,meta$FRiP,method = 'pearson')
corr_FRiP_4 = cor(umaps$X4,meta$FRiP,method = 'pearson')
corr_FRiP_5 = cor(umaps$X5,meta$FRiP,method = 'pearson')
corr_FRiP_6 = cor(umaps$X6,meta$FRiP,method = 'pearson')

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
group_list2 = meta$unique_reads
group_list3 = meta$Peaks
group_list4 = meta$NSC
group_list5 = meta$RSC
group_list6 = meta$FRiP
group_list7 = meta$sequencing_depth
group_list8 = meta$batch

colData=data.frame(row.names = colnames(group_CA3),
                   condition=group_list1,
                   reads=group_list2,
                   peaks=group_list3,
                   NSC=group_list4,
                   RSC = group_list5,
                   FRiP = group_list6,
                   sequencing_depth = group_list7,
                   batch = group_list8)

colData$condition = factor(colData$condition,levels = c('C','A'))
colData$reads =as.numeric(colData$reads)
colData$peaks=as.numeric(colData$peaks)
colData$NSC=as.numeric(colData$NSC)
colData$RSC=as.numeric(colData$RSC)
colData$FRiP=as.numeric(colData$FRiP)
colData$sequencing_depth=as.numeric(colData$sequencing_depth)
colData$batch=as.factor(colData$batch)

dds=DESeqDataSetFromMatrix(countData = group_CA3,
                           colData = colData,
                           design = ~ condition + NSC)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","C","A"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) >0.6),]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_CA3[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0 & nrDEG_glia$log2FoldChange > 0.5),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0 & nrDEG_glia$log2FoldChange < -0.5),]
diff_peak = peak_info[rownames(DEG_matrix),]
new_DEG_matrix = bind_cols(diff_peak,DEG_matrix)
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/CA3")
all_glia = group_CA3[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_CA3.csv')
write.csv(new_DEG_matrix,'Differential_CA3.csv')
# ###############################B vs. D_1###############################################################
meta = bind_rows(factors[29:40,1:8],factors[53:64,1:8])
meta$condition_num = c(rep(1,12),rep(2,12))
group_BD1 = data.frame(def[,29:40],def[,53:64])
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
cor_con_peaks = cor(meta$condition_num,meta$Peaks,method = 'pearson')
cor_con_NSC = cor(meta$condition_num,meta$NSC,method = 'pearson')
cor_con_RSC = cor(meta$condition_num,meta$RSC,method = 'pearson')
cor_con_FRiP = cor(meta$condition_num,meta$FRiP,method = 'pearson')
cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
cor_con_uni_reads = cor(meta$condition_num,meta$unique_reads,method = 'pearson')
# FRiP
# FRiP NSC -UMAP
corr_uni_reads_1 = cor(umaps$X1,meta$unique_reads,method = 'pearson')
corr_uni_reads_2 = cor(umaps$X2,meta$unique_reads,method = 'pearson')
corr_uni_reads_3 = cor(umaps$X3,meta$unique_reads,method = 'pearson')
corr_uni_reads_4 = cor(umaps$X4,meta$unique_reads,method = 'pearson')
corr_uni_reads_5 = cor(umaps$X5,meta$unique_reads,method = 'pearson')
corr_uni_reads_6 = cor(umaps$X6,meta$unique_reads,method = 'pearson')

corr_peaks_1 = cor(umaps$X1,meta$Peaks,method = 'pearson')
corr_peaks_2 = cor(umaps$X2,meta$Peaks,method = 'pearson')
corr_peaks_3 = cor(umaps$X3,meta$Peaks,method = 'pearson')
corr_peaks_4 = cor(umaps$X4,meta$Peaks,method = 'pearson')
corr_peaks_5 = cor(umaps$X5,meta$Peaks,method = 'pearson')
corr_peaks_6 = cor(umaps$X6,meta$Peaks,method = 'pearson')

corr_NSC_1 = cor(umaps$X1,meta$NSC,method = 'pearson')
corr_NSC_2 = cor(umaps$X2,meta$NSC,method = 'pearson')
corr_NSC_3 = cor(umaps$X3,meta$NSC,method = 'pearson')
corr_NSC_4 = cor(umaps$X4,meta$NSC,method = 'pearson')
corr_NSC_5 = cor(umaps$X5,meta$NSC,method = 'pearson')
corr_NSC_6 = cor(umaps$X6,meta$NSC,method = 'pearson')

corr_RSC_1 = cor(umaps$X1,meta$RSC,method = 'pearson')
corr_RSC_2 = cor(umaps$X2,meta$RSC,method = 'pearson')
corr_RSC_3 = cor(umaps$X3,meta$RSC,method = 'pearson')
corr_RSC_4 = cor(umaps$X4,meta$RSC,method = 'pearson')
corr_RSC_5 = cor(umaps$X5,meta$RSC,method = 'pearson')
corr_RSC_6 = cor(umaps$X6,meta$RSC,method = 'pearson')

corr_FRiP_1 = cor(umaps$X1,meta$FRiP,method = 'pearson')
corr_FRiP_2 = cor(umaps$X2,meta$FRiP,method = 'pearson')
corr_FRiP_3 = cor(umaps$X3,meta$FRiP,method = 'pearson')
corr_FRiP_4 = cor(umaps$X4,meta$FRiP,method = 'pearson')
corr_FRiP_5 = cor(umaps$X5,meta$FRiP,method = 'pearson')
corr_FRiP_6 = cor(umaps$X6,meta$FRiP,method = 'pearson')

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
group_list2 = meta$unique_reads
group_list3 = meta$Peaks
group_list4 = meta$NSC
group_list5 = meta$RSC
group_list6 = meta$FRiP
group_list7 = meta$sequencing_depth
group_list8 = meta$batch

colData=data.frame(row.names = colnames(group_BD1),
                   condition=group_list1,
                   reads=group_list2,
                   peaks=group_list3,
                   NSC=group_list4,
                   RSC = group_list5,
                   FRiP = group_list6,
                   sequencing_depth = group_list7,
                   batch = group_list8)

colData$condition = factor(colData$condition,levels = c('B','D'))
colData$reads =as.numeric(colData$reads)
colData$peaks=as.numeric(colData$peaks)
colData$NSC=as.numeric(colData$NSC)
colData$RSC=as.numeric(colData$RSC)
colData$FRiP=as.numeric(colData$FRiP)
colData$sequencing_depth=as.numeric(colData$sequencing_depth)
colData$batch=as.factor(colData$batch)

dds=DESeqDataSetFromMatrix(countData = group_BD1,
                           colData = colData,
                           design = ~ condition + NSC + FRiP + RSC)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","B","D"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) >0.6),]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_BD1[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0 & nrDEG_glia$log2FoldChange > 0.5),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0 & nrDEG_glia$log2FoldChange < -0.5),]
diff_peak = peak_info[rownames(DEG_matrix),]
new_DEG_matrix = bind_cols(diff_peak,DEG_matrix)
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/BD1")
all_glia = group_BD1[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_BD1.csv')
write.csv(new_DEG_matrix,'Differential_BD1.csv')
# ###############################B vs. D_2###############################################################
meta = bind_rows(factors[29:40,1:8],factors[65:70,1:8])
meta$condition_num = c(rep(1,12),rep(2,6))
group_BD2 = data.frame(def[,29:40],def[,65:70])
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
cor_con_peaks = cor(meta$condition_num,meta$Peaks,method = 'pearson')
cor_con_NSC = cor(meta$condition_num,meta$NSC,method = 'pearson')
cor_con_RSC = cor(meta$condition_num,meta$RSC,method = 'pearson')
cor_con_FRiP = cor(meta$condition_num,meta$FRiP,method = 'pearson')
cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
cor_con_uni_reads = cor(meta$condition_num,meta$unique_reads,method = 'pearson')
# peaks RSC
# peaks RSC reads -UMAP
corr_uni_reads_1 = cor(umaps$X1,meta$unique_reads,method = 'pearson')
corr_uni_reads_2 = cor(umaps$X2,meta$unique_reads,method = 'pearson')
corr_uni_reads_3 = cor(umaps$X3,meta$unique_reads,method = 'pearson')
corr_uni_reads_4 = cor(umaps$X4,meta$unique_reads,method = 'pearson')
corr_uni_reads_5 = cor(umaps$X5,meta$unique_reads,method = 'pearson')
corr_uni_reads_6 = cor(umaps$X6,meta$unique_reads,method = 'pearson')

corr_peaks_1 = cor(umaps$X1,meta$Peaks,method = 'pearson')
corr_peaks_2 = cor(umaps$X2,meta$Peaks,method = 'pearson')
corr_peaks_3 = cor(umaps$X3,meta$Peaks,method = 'pearson')
corr_peaks_4 = cor(umaps$X4,meta$Peaks,method = 'pearson')
corr_peaks_5 = cor(umaps$X5,meta$Peaks,method = 'pearson')
corr_peaks_6 = cor(umaps$X6,meta$Peaks,method = 'pearson')

corr_NSC_1 = cor(umaps$X1,meta$NSC,method = 'pearson')
corr_NSC_2 = cor(umaps$X2,meta$NSC,method = 'pearson')
corr_NSC_3 = cor(umaps$X3,meta$NSC,method = 'pearson')
corr_NSC_4 = cor(umaps$X4,meta$NSC,method = 'pearson')
corr_NSC_5 = cor(umaps$X5,meta$NSC,method = 'pearson')
corr_NSC_6 = cor(umaps$X6,meta$NSC,method = 'pearson')

corr_RSC_1 = cor(umaps$X1,meta$RSC,method = 'pearson')
corr_RSC_2 = cor(umaps$X2,meta$RSC,method = 'pearson')
corr_RSC_3 = cor(umaps$X3,meta$RSC,method = 'pearson')
corr_RSC_4 = cor(umaps$X4,meta$RSC,method = 'pearson')
corr_RSC_5 = cor(umaps$X5,meta$RSC,method = 'pearson')
corr_RSC_6 = cor(umaps$X6,meta$RSC,method = 'pearson')

corr_FRiP_1 = cor(umaps$X1,meta$FRiP,method = 'pearson')
corr_FRiP_2 = cor(umaps$X2,meta$FRiP,method = 'pearson')
corr_FRiP_3 = cor(umaps$X3,meta$FRiP,method = 'pearson')
corr_FRiP_4 = cor(umaps$X4,meta$FRiP,method = 'pearson')
corr_FRiP_5 = cor(umaps$X5,meta$FRiP,method = 'pearson')
corr_FRiP_6 = cor(umaps$X6,meta$FRiP,method = 'pearson')

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
group_list2 = meta$unique_reads
group_list3 = meta$Peaks
group_list4 = meta$NSC
group_list5 = meta$RSC
group_list6 = meta$FRiP
group_list7 = meta$sequencing_depth
group_list8 = meta$batch

colData=data.frame(row.names = colnames(group_BD2),
                   condition=group_list1,
                   reads=group_list2,
                   peaks=group_list3,
                   NSC=group_list4,
                   RSC = group_list5,
                   FRiP = group_list6,
                   sequencing_depth = group_list7,
                   batch = group_list8)

colData$condition = factor(colData$condition,levels = c('B','D'))
colData$reads =as.numeric(colData$reads)
colData$peaks=as.numeric(colData$peaks)
colData$NSC=as.numeric(colData$NSC)
colData$RSC=as.numeric(colData$RSC)
colData$FRiP=as.numeric(colData$FRiP)
colData$sequencing_depth=as.numeric(colData$sequencing_depth)
colData$batch=as.factor(colData$batch)

dds=DESeqDataSetFromMatrix(countData = group_BD2,
                           colData = colData,
                           design = ~ condition + reads  + NSC + peaks)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","B","D"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) >0.6),]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_BD2[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0 & nrDEG_glia$log2FoldChange > 0.5),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0 & nrDEG_glia$log2FoldChange < -0.5),]
diff_peak = peak_info[rownames(DEG_matrix),]
new_DEG_matrix = bind_cols(diff_peak,DEG_matrix)
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/BD2")
all_glia = group_BD2[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_BD2.csv')
write.csv(new_DEG_matrix,'Differential_BD2.csv')
# ###############################B vs. D_3###############################################################
meta = bind_rows(factors[29:40,1:8],factors[71:82,1:8])
meta$condition_num = c(rep(1,12),rep(2,12))
group_BD3 = data.frame(def[,29:40],def[,71:82])
miss_1 <- c()
for (i in 1:nrow(group_BD3)) {
  if(length(which(group_BD3[i,1:12]<10)) > 0.5*ncol(group_BD3[,1:12])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(group_BD3)) {
  if(length(which(group_BD3[i,13:24]<10)) > 0.5*ncol(group_BD3[,13:24])) miss_2 <- append(miss_2,i)
}
miss = union(miss_1,miss_2)
group_BD3 = group_BD3[-miss,]
t_group_BD3 = t(group_BD3)
umap_results <- umap(t_group_BD3,n_neighbors=11,min_dist=0.1,n_components = 6)
umaps = data.frame(umap_results[["layout"]])
umaps$condition = as.factor(meta$condition)
ggplot(umaps,mapping = aes(x=X1,y=X2,col=condition)) +
  geom_point(aes(colour = condition))

cor_con_batch = cor(meta$condition_num,meta$batch,method = 'pearson')
cor_con_peaks = cor(meta$condition_num,meta$Peaks,method = 'pearson')
cor_con_NSC = cor(meta$condition_num,meta$NSC,method = 'pearson')
cor_con_RSC = cor(meta$condition_num,meta$RSC,method = 'pearson')
cor_con_FRiP = cor(meta$condition_num,meta$FRiP,method = 'pearson')
cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
cor_con_uni_reads = cor(meta$condition_num,meta$unique_reads,method = 'pearson')
# NSC peaks reads
# NSC peaks reads -UMAP
corr_uni_reads_1 = cor(umaps$X1,meta$unique_reads,method = 'pearson')
corr_uni_reads_2 = cor(umaps$X2,meta$unique_reads,method = 'pearson')
corr_uni_reads_3 = cor(umaps$X3,meta$unique_reads,method = 'pearson')
corr_uni_reads_4 = cor(umaps$X4,meta$unique_reads,method = 'pearson')
corr_uni_reads_5 = cor(umaps$X5,meta$unique_reads,method = 'pearson')
corr_uni_reads_6 = cor(umaps$X6,meta$unique_reads,method = 'pearson')

corr_peaks_1 = cor(umaps$X1,meta$Peaks,method = 'pearson')
corr_peaks_2 = cor(umaps$X2,meta$Peaks,method = 'pearson')
corr_peaks_3 = cor(umaps$X3,meta$Peaks,method = 'pearson')
corr_peaks_4 = cor(umaps$X4,meta$Peaks,method = 'pearson')
corr_peaks_5 = cor(umaps$X5,meta$Peaks,method = 'pearson')
corr_peaks_6 = cor(umaps$X6,meta$Peaks,method = 'pearson')

corr_NSC_1 = cor(umaps$X1,meta$NSC,method = 'pearson')
corr_NSC_2 = cor(umaps$X2,meta$NSC,method = 'pearson')
corr_NSC_3 = cor(umaps$X3,meta$NSC,method = 'pearson')
corr_NSC_4 = cor(umaps$X4,meta$NSC,method = 'pearson')
corr_NSC_5 = cor(umaps$X5,meta$NSC,method = 'pearson')
corr_NSC_6 = cor(umaps$X6,meta$NSC,method = 'pearson')

corr_RSC_1 = cor(umaps$X1,meta$RSC,method = 'pearson')
corr_RSC_2 = cor(umaps$X2,meta$RSC,method = 'pearson')
corr_RSC_3 = cor(umaps$X3,meta$RSC,method = 'pearson')
corr_RSC_4 = cor(umaps$X4,meta$RSC,method = 'pearson')
corr_RSC_5 = cor(umaps$X5,meta$RSC,method = 'pearson')
corr_RSC_6 = cor(umaps$X6,meta$RSC,method = 'pearson')

corr_FRiP_1 = cor(umaps$X1,meta$FRiP,method = 'pearson')
corr_FRiP_2 = cor(umaps$X2,meta$FRiP,method = 'pearson')
corr_FRiP_3 = cor(umaps$X3,meta$FRiP,method = 'pearson')
corr_FRiP_4 = cor(umaps$X4,meta$FRiP,method = 'pearson')
corr_FRiP_5 = cor(umaps$X5,meta$FRiP,method = 'pearson')
corr_FRiP_6 = cor(umaps$X6,meta$FRiP,method = 'pearson')

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
group_list2 = meta$unique_reads
group_list3 = meta$Peaks
group_list4 = meta$NSC
group_list5 = meta$RSC
group_list6 = meta$FRiP
group_list7 = meta$sequencing_depth
group_list8 = meta$batch

colData=data.frame(row.names = colnames(group_BD3),
                   condition=group_list1,
                   reads=group_list2,
                   peaks=group_list3,
                   NSC=group_list4,
                   RSC = group_list5,
                   FRiP = group_list6,
                   sequencing_depth = group_list7,
                   batch = group_list8)

colData$condition = factor(colData$condition,levels = c('B','D'))
colData$reads =as.numeric(colData$reads)
colData$peaks=as.numeric(colData$peaks)
colData$NSC=as.numeric(colData$NSC)
colData$RSC=as.numeric(colData$RSC)
colData$FRiP=as.numeric(colData$FRiP)
colData$sequencing_depth=as.numeric(colData$sequencing_depth)
colData$batch=as.factor(colData$batch)

dds=DESeqDataSetFromMatrix(countData = group_BD3,
                           colData = colData,
                           design = ~ condition + NSC + reads)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","B","D"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) >0.6),]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_BD3[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0 & nrDEG_glia$log2FoldChange > 0.5),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0 & nrDEG_glia$log2FoldChange < -0.5),]
diff_peak = peak_info[rownames(DEG_matrix),]
new_DEG_matrix = bind_cols(diff_peak,DEG_matrix)
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/BD3")
all_glia = group_BD3[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_BD3.csv')
write.csv(new_DEG_matrix,'Differential_BD3.csv')
# #########################EFG######################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/count_files/")
raw_def = read.csv('promoter_efg_count.csv')
rownames(raw_def) = raw_def[,4]
def = raw_def[,5:ncol(raw_def)]
peak_info = raw_def[,1:4]
# ###############################F vs. G###############################################################
meta = bind_rows(factors[95:106,1:8],factors[107:118,1:8])
meta$condition_num = c(rep(1,12),rep(2,12))
group_FG = data.frame(def[,13:24],def[,25:36])
miss_1 <- c()
for (i in 1:nrow(group_FG)) {
  if(length(which(group_FG[i,1:12]<10)) > 0.5*ncol(group_FG[,1:12])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(group_FG)) {
  if(length(which(group_FG[i,12:24]<10)) > 0.5*ncol(group_FG[,12:24])) miss_2 <- append(miss_2,i)
}
miss = union(miss_1,miss_2)
group_FG = group_FG[-miss,]
t_group_FG = t(group_FG)
umap_results <- umap(t_group_FG,n_neighbors=12,min_dist=0.1,n_components = 6)
umaps = data.frame(umap_results[["layout"]])
umaps$condition = as.factor(meta$condition)
ggplot(umaps,mapping = aes(x=X1,y=X2,col=condition)) +
  geom_point(aes(colour = condition))

cor_con_batch = cor(meta$condition_num,meta$batch,method = 'pearson')
cor_con_peaks = cor(meta$condition_num,meta$Peaks,method = 'pearson')
cor_con_NSC = cor(meta$condition_num,meta$NSC,method = 'pearson')
cor_con_RSC = cor(meta$condition_num,meta$RSC,method = 'pearson')
cor_con_FRiP = cor(meta$condition_num,meta$FRiP,method = 'pearson')
cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
cor_con_uni_reads = cor(meta$condition_num,meta$unique_reads,method = 'pearson')
# peaks
# peaks reads -UMAP
corr_uni_reads_1 = cor(umaps$X1,meta$unique_reads,method = 'pearson')
corr_uni_reads_2 = cor(umaps$X2,meta$unique_reads,method = 'pearson')
corr_uni_reads_3 = cor(umaps$X3,meta$unique_reads,method = 'pearson')
corr_uni_reads_4 = cor(umaps$X4,meta$unique_reads,method = 'pearson')
corr_uni_reads_5 = cor(umaps$X5,meta$unique_reads,method = 'pearson')
corr_uni_reads_6 = cor(umaps$X6,meta$unique_reads,method = 'pearson')

corr_peaks_1 = cor(umaps$X1,meta$Peaks,method = 'pearson')
corr_peaks_2 = cor(umaps$X2,meta$Peaks,method = 'pearson')
corr_peaks_3 = cor(umaps$X3,meta$Peaks,method = 'pearson')
corr_peaks_4 = cor(umaps$X4,meta$Peaks,method = 'pearson')
corr_peaks_5 = cor(umaps$X5,meta$Peaks,method = 'pearson')
corr_peaks_6 = cor(umaps$X6,meta$Peaks,method = 'pearson')

corr_NSC_1 = cor(umaps$X1,meta$NSC,method = 'pearson')
corr_NSC_2 = cor(umaps$X2,meta$NSC,method = 'pearson')
corr_NSC_3 = cor(umaps$X3,meta$NSC,method = 'pearson')
corr_NSC_4 = cor(umaps$X4,meta$NSC,method = 'pearson')
corr_NSC_5 = cor(umaps$X5,meta$NSC,method = 'pearson')
corr_NSC_6 = cor(umaps$X6,meta$NSC,method = 'pearson')

corr_RSC_1 = cor(umaps$X1,meta$RSC,method = 'pearson')
corr_RSC_2 = cor(umaps$X2,meta$RSC,method = 'pearson')
corr_RSC_3 = cor(umaps$X3,meta$RSC,method = 'pearson')
corr_RSC_4 = cor(umaps$X4,meta$RSC,method = 'pearson')
corr_RSC_5 = cor(umaps$X5,meta$RSC,method = 'pearson')
corr_RSC_6 = cor(umaps$X6,meta$RSC,method = 'pearson')

corr_FRiP_1 = cor(umaps$X1,meta$FRiP,method = 'pearson')
corr_FRiP_2 = cor(umaps$X2,meta$FRiP,method = 'pearson')
corr_FRiP_3 = cor(umaps$X3,meta$FRiP,method = 'pearson')
corr_FRiP_4 = cor(umaps$X4,meta$FRiP,method = 'pearson')
corr_FRiP_5 = cor(umaps$X5,meta$FRiP,method = 'pearson')
corr_FRiP_6 = cor(umaps$X6,meta$FRiP,method = 'pearson')

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
group_list2 = meta$unique_reads
group_list3 = meta$Peaks
group_list4 = meta$NSC
group_list5 = meta$RSC
group_list6 = meta$FRiP
group_list7 = meta$sequencing_depth
group_list8 = meta$batch

colData=data.frame(row.names = colnames(group_FG),
                   condition=group_list1,
                   reads=group_list2,
                   peaks=group_list3,
                   NSC=group_list4,
                   RSC = group_list5,
                   FRiP = group_list6,
                   sequencing_depth = group_list7)

colData$condition = factor(colData$condition,levels = c('F','G'))
colData$reads =as.numeric(colData$reads)
colData$peaks=as.numeric(colData$peaks)
colData$NSC=as.numeric(colData$NSC)
colData$RSC=as.numeric(colData$RSC)
colData$FRiP=as.numeric(colData$FRiP)
colData$sequencing_depth=as.numeric(colData$sequencing_depth)

dds=DESeqDataSetFromMatrix(countData = group_FG,
                           colData = colData,
                           design = ~ condition + peaks)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","F","G"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.6),]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_FG[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0 & nrDEG_glia$log2FoldChange > 0.5),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0 & nrDEG_glia$log2FoldChange < -0.5),]
diff_peak = peak_info[rownames(DEG_matrix),]
new_DEG_matrix = bind_cols(diff_peak,DEG_matrix)
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/FG")
all_glia = group_FG[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_FG.csv')
write.csv(new_DEG_matrix,'Differential_FG.csv')
# ###############################E vs. G###############################################################
meta = bind_rows(factors[83:94,1:8],factors[107:118,1:8])
meta$condition_num = c(rep(1,12),rep(2,12))
group_EG = data.frame(def[,1:12],def[,25:36])
miss_1 <- c()
for (i in 1:nrow(group_EG)) {
  if(length(which(group_EG[i,1:12]<10)) > 0.5*ncol(group_EG[,1:12])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(group_EG)) {
  if(length(which(group_EG[i,13:24]<10)) > 0.5*ncol(group_EG[,13:24])) miss_2 <- append(miss_2,i)
}
miss = union(miss_1,miss_2)
group_EG = group_EG[-miss,]
t_group_EG = t(group_EG)
umap_results <- umap(t_group_EG,n_neighbors=12,min_dist=0.1,n_components = 6)
umaps = data.frame(umap_results[["layout"]])
umaps$condition = as.factor(meta$condition)
ggplot(umaps,mapping = aes(x=X1,y=X2,col=condition)) +
  geom_point(aes(colour = condition))

cor_con_batch = cor(meta$condition_num,meta$batch,method = 'pearson')
cor_con_peaks = cor(meta$condition_num,meta$Peaks,method = 'pearson')
cor_con_NSC = cor(meta$condition_num,meta$NSC,method = 'pearson')
cor_con_RSC = cor(meta$condition_num,meta$RSC,method = 'pearson')
cor_con_FRiP = cor(meta$condition_num,meta$FRiP,method = 'pearson')
cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
cor_con_uni_reads = cor(meta$condition_num,meta$unique_reads,method = 'pearson')
# peaks
# NSC -UMAP
corr_uni_reads_1 = cor(umaps$X1,meta$unique_reads,method = 'pearson')
corr_uni_reads_2 = cor(umaps$X2,meta$unique_reads,method = 'pearson')
corr_uni_reads_3 = cor(umaps$X3,meta$unique_reads,method = 'pearson')
corr_uni_reads_4 = cor(umaps$X4,meta$unique_reads,method = 'pearson')
corr_uni_reads_5 = cor(umaps$X5,meta$unique_reads,method = 'pearson')
corr_uni_reads_6 = cor(umaps$X6,meta$unique_reads,method = 'pearson')

corr_peaks_1 = cor(umaps$X1,meta$Peaks,method = 'pearson')
corr_peaks_2 = cor(umaps$X2,meta$Peaks,method = 'pearson')
corr_peaks_3 = cor(umaps$X3,meta$Peaks,method = 'pearson')
corr_peaks_4 = cor(umaps$X4,meta$Peaks,method = 'pearson')
corr_peaks_5 = cor(umaps$X5,meta$Peaks,method = 'pearson')
corr_peaks_6 = cor(umaps$X6,meta$Peaks,method = 'pearson')

corr_NSC_1 = cor(umaps$X1,meta$NSC,method = 'pearson')
corr_NSC_2 = cor(umaps$X2,meta$NSC,method = 'pearson')
corr_NSC_3 = cor(umaps$X3,meta$NSC,method = 'pearson')
corr_NSC_4 = cor(umaps$X4,meta$NSC,method = 'pearson')
corr_NSC_5 = cor(umaps$X5,meta$NSC,method = 'pearson')
corr_NSC_6 = cor(umaps$X6,meta$NSC,method = 'pearson')

corr_RSC_1 = cor(umaps$X1,meta$RSC,method = 'pearson')
corr_RSC_2 = cor(umaps$X2,meta$RSC,method = 'pearson')
corr_RSC_3 = cor(umaps$X3,meta$RSC,method = 'pearson')
corr_RSC_4 = cor(umaps$X4,meta$RSC,method = 'pearson')
corr_RSC_5 = cor(umaps$X5,meta$RSC,method = 'pearson')
corr_RSC_6 = cor(umaps$X6,meta$RSC,method = 'pearson')

corr_FRiP_1 = cor(umaps$X1,meta$FRiP,method = 'pearson')
corr_FRiP_2 = cor(umaps$X2,meta$FRiP,method = 'pearson')
corr_FRiP_3 = cor(umaps$X3,meta$FRiP,method = 'pearson')
corr_FRiP_4 = cor(umaps$X4,meta$FRiP,method = 'pearson')
corr_FRiP_5 = cor(umaps$X5,meta$FRiP,method = 'pearson')
corr_FRiP_6 = cor(umaps$X6,meta$FRiP,method = 'pearson')

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
group_list2 = meta$unique_reads
group_list3 = meta$Peaks
group_list4 = meta$NSC
group_list5 = meta$RSC
group_list6 = meta$FRiP
group_list7 = meta$sequencing_depth
group_list8 = meta$batch

colData=data.frame(row.names = colnames(group_EG),
                   condition=group_list1,
                   reads=group_list2,
                   peaks=group_list3,
                   NSC=group_list4,
                   RSC = group_list5,
                   FRiP = group_list6,
                   sequencing_depth = group_list7,
                   batch = group_list8)

colData$condition = factor(colData$condition,levels = c('E','G'))
colData$reads =as.numeric(colData$reads)
colData$peaks=as.numeric(colData$peaks)
colData$NSC=as.numeric(colData$NSC)
colData$RSC=as.numeric(colData$RSC)
colData$FRiP=as.numeric(colData$FRiP)
colData$sequencing_depth=as.numeric(colData$sequencing_depth)
colData$batch=as.factor(colData$batch)

dds=DESeqDataSetFromMatrix(countData = group_EG,
                           colData = colData,
                           design = ~ condition + peaks)

# dds=DESeqDataSetFromMatrix(countData = group_EG,
#                            colData = colData,
#                            design = ~ condition)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","E","G"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.6),]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_EG[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0 & nrDEG_glia$log2FoldChange > 0.5),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0 & nrDEG_glia$log2FoldChange < -0.5),]
diff_peak = peak_info[rownames(DEG_matrix),]
new_DEG_matrix = bind_cols(diff_peak,DEG_matrix)
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/EG")
all_glia = group_EG[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_EG.csv')
write.csv(new_DEG_matrix,'Differential_EG.csv')
###############################F vs. E###############################################################
meta = bind_rows(factors[95:106,1:8],factors[83:94,1:8])
meta$condition_num = c(rep(1,12),rep(2,12))
group_FE = data.frame(def[,13:24],def[,1:12])
miss_1 <- c()
for (i in 1:nrow(group_FE)) {
  if(length(which(group_FE[i,1:12]<10)) > 0.5*ncol(group_FE[,1:12])) miss_1 <- append(miss_1,i)
}
miss_2 <- c()
for (i in 1:nrow(group_FE)) {
  if(length(which(group_FE[i,13:24]<10)) > 0.5*ncol(group_FE[,13:24])) miss_2 <- append(miss_2,i)
}
miss = union(miss_1,miss_2)
group_FE = group_FE[-miss,]
t_group_FE = t(group_FE)
umap_results <- umap(t_group_FE,n_neighbors=12,min_dist=0.1,n_components = 6)
umaps = data.frame(umap_results[["layout"]])
umaps$condition = as.factor(meta$condition)
ggplot(umaps,mapping = aes(x=X1,y=X2,col=condition)) +
  geom_point(aes(colour = condition))

cor_con_batch = cor(meta$condition_num,meta$batch,method = 'pearson')
cor_con_peaks = cor(meta$condition_num,meta$Peaks,method = 'pearson')
cor_con_NSC = cor(meta$condition_num,meta$NSC,method = 'pearson')
cor_con_RSC = cor(meta$condition_num,meta$RSC,method = 'pearson')
cor_con_FRiP = cor(meta$condition_num,meta$FRiP,method = 'pearson')
cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
cor_con_uni_reads = cor(meta$condition_num,meta$unique_reads,method = 'pearson')
#
#  NSC -UMAP
corr_uni_reads_1 = cor(umaps$X1,meta$unique_reads,method = 'pearson')
corr_uni_reads_2 = cor(umaps$X2,meta$unique_reads,method = 'pearson')
corr_uni_reads_3 = cor(umaps$X3,meta$unique_reads,method = 'pearson')
corr_uni_reads_4 = cor(umaps$X4,meta$unique_reads,method = 'pearson')
corr_uni_reads_5 = cor(umaps$X5,meta$unique_reads,method = 'pearson')
corr_uni_reads_6 = cor(umaps$X6,meta$unique_reads,method = 'pearson')

corr_peaks_1 = cor(umaps$X1,meta$Peaks,method = 'pearson')
corr_peaks_2 = cor(umaps$X2,meta$Peaks,method = 'pearson')
corr_peaks_3 = cor(umaps$X3,meta$Peaks,method = 'pearson')
corr_peaks_4 = cor(umaps$X4,meta$Peaks,method = 'pearson')
corr_peaks_5 = cor(umaps$X5,meta$Peaks,method = 'pearson')
corr_peaks_6 = cor(umaps$X6,meta$Peaks,method = 'pearson')

corr_NSC_1 = cor(umaps$X1,meta$NSC,method = 'pearson')
corr_NSC_2 = cor(umaps$X2,meta$NSC,method = 'pearson')
corr_NSC_3 = cor(umaps$X3,meta$NSC,method = 'pearson')
corr_NSC_4 = cor(umaps$X4,meta$NSC,method = 'pearson')
corr_NSC_5 = cor(umaps$X5,meta$NSC,method = 'pearson')
corr_NSC_6 = cor(umaps$X6,meta$NSC,method = 'pearson')

corr_RSC_1 = cor(umaps$X1,meta$RSC,method = 'pearson')
corr_RSC_2 = cor(umaps$X2,meta$RSC,method = 'pearson')
corr_RSC_3 = cor(umaps$X3,meta$RSC,method = 'pearson')
corr_RSC_4 = cor(umaps$X4,meta$RSC,method = 'pearson')
corr_RSC_5 = cor(umaps$X5,meta$RSC,method = 'pearson')
corr_RSC_6 = cor(umaps$X6,meta$RSC,method = 'pearson')

corr_FRiP_1 = cor(umaps$X1,meta$FRiP,method = 'pearson')
corr_FRiP_2 = cor(umaps$X2,meta$FRiP,method = 'pearson')
corr_FRiP_3 = cor(umaps$X3,meta$FRiP,method = 'pearson')
corr_FRiP_4 = cor(umaps$X4,meta$FRiP,method = 'pearson')
corr_FRiP_5 = cor(umaps$X5,meta$FRiP,method = 'pearson')
corr_FRiP_6 = cor(umaps$X6,meta$FRiP,method = 'pearson')

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
group_list2 = meta$unique_reads
group_list3 = meta$Peaks
group_list4 = meta$NSC
group_list5 = meta$RSC
group_list6 = meta$FRiP
group_list7 = meta$sequencing_depth
group_list8 = meta$batch

colData=data.frame(row.names = colnames(group_FE),
                   condition=group_list1,
                   reads=group_list2,
                   peaks=group_list3,
                   NSC=group_list4,
                   RSC = group_list5,
                   FRiP = group_list6,
                   sequencing_depth = group_list7,
                   batch = group_list8)

colData$condition = factor(colData$condition,levels = c('F','E'))
colData$reads =as.numeric(colData$reads)
colData$peaks=as.numeric(colData$peaks)
colData$NSC=as.numeric(colData$NSC)
colData$RSC=as.numeric(colData$RSC)
colData$FRiP=as.numeric(colData$FRiP)
colData$sequencing_depth=as.numeric(colData$sequencing_depth)
colData$batch=as.factor(colData$batch)

# dds=DESeqDataSetFromMatrix(countData = group_FE,
#                            colData = colData,
#                            design = ~ condition)

dds=DESeqDataSetFromMatrix(countData = group_FE,
                           colData = colData,
                           design = ~ condition + NSC + FRiP)
#tmp = model.matrix(~group_list + gender + Age + PMI, colData)

dds=DESeq(dds)

res_glia=results(dds,
                 contrast = c("condition","F","E"))
resOrdered_glia=res_glia[order(res_glia$padj),]
DEG_glia=as.data.frame(resOrdered_glia)
DEG_glia=na.omit(DEG_glia)
nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.6),]
choose_gene_neuron=rownames(nrDEG_glia)
choose_matrix=group_FE[choose_gene_neuron,]
DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0 & nrDEG_glia$log2FoldChange > 0.5),]
down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0 & nrDEG_glia$log2FoldChange < -0.5),]
diff_peak = peak_info[rownames(DEG_matrix),]
new_DEG_matrix = bind_cols(diff_peak,DEG_matrix)
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/FE")
all_glia = group_FE[rownames(DEG_glia),]
all_matrix = data.frame(DEG_glia,all_glia)
write.csv(all_matrix,'All_FE.csv')
write.csv(new_DEG_matrix,'Differential_FE.csv')

# ###############################D_1 vs. A_1###############################################################
# meta = bind_rows(factors[54:65,1:8],factors[1:11,1:8])
# meta$condition_num = c(rep(1,12),rep(2,11))
# group_D1A1 = data.frame(def[,54:65],def[,1:11])
# miss_1 <- c()
# for (i in 1:nrow(group_D1A1)) {
#   if(length(which(group_D1A1[i,1:12]<10)) > 0.5*ncol(group_D1A1[,1:12])) miss_1 <- append(miss_1,i)
# }
# miss_2 <- c()
# for (i in 1:nrow(group_D1A1)) {
#   if(length(which(group_D1A1[i,13:23]<10)) > 0.5*ncol(group_D1A1[,13:23])) miss_2 <- append(miss_2,i)
# }
# miss = union(miss_1,miss_2)
# group_D1A1 = group_D1A1[-miss,]
# t_group_D1A1 = t(group_D1A1)
# umap_results <- umap(t_group_D1A1,n_neighbors=11,min_dist=0.1,n_components = 6)
# umaps = data.frame(umap_results[["layout"]])
# umaps$condition = as.factor(meta$condition)
# ggplot(umaps,mapping = aes(x=X1,y=X2,col=condition)) +
#   geom_point(aes(colour = condition))
# 
# cor_con_batch = cor(meta$condition_num,meta$batch,method = 'pearson')
# cor_con_peaks = cor(meta$condition_num,meta$Peaks,method = 'pearson')
# cor_con_NSC = cor(meta$condition_num,meta$NSC,method = 'pearson')
# cor_con_RSC = cor(meta$condition_num,meta$RSC,method = 'pearson')
# cor_con_FRiP = cor(meta$condition_num,meta$FRiP,method = 'pearson')
# cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
# cor_con_uni_reads = cor(meta$condition_num,meta$unique_reads,method = 'pearson')
# #  
# # NSC RSC -UMAP
# corr_uni_reads_1 = cor(umaps$X1,meta$unique_reads,method = 'pearson')
# corr_uni_reads_2 = cor(umaps$X2,meta$unique_reads,method = 'pearson')
# corr_uni_reads_3 = cor(umaps$X3,meta$unique_reads,method = 'pearson')
# corr_uni_reads_4 = cor(umaps$X4,meta$unique_reads,method = 'pearson')
# corr_uni_reads_5 = cor(umaps$X5,meta$unique_reads,method = 'pearson')
# corr_uni_reads_6 = cor(umaps$X6,meta$unique_reads,method = 'pearson')
# 
# corr_peaks_1 = cor(umaps$X1,meta$Peaks,method = 'pearson')
# corr_peaks_2 = cor(umaps$X2,meta$Peaks,method = 'pearson')
# corr_peaks_3 = cor(umaps$X3,meta$Peaks,method = 'pearson')
# corr_peaks_4 = cor(umaps$X4,meta$Peaks,method = 'pearson')
# corr_peaks_5 = cor(umaps$X5,meta$Peaks,method = 'pearson')
# corr_peaks_6 = cor(umaps$X6,meta$Peaks,method = 'pearson')
# 
# corr_NSC_1 = cor(umaps$X1,meta$NSC,method = 'pearson')
# corr_NSC_2 = cor(umaps$X2,meta$NSC,method = 'pearson')
# corr_NSC_3 = cor(umaps$X3,meta$NSC,method = 'pearson')
# corr_NSC_4 = cor(umaps$X4,meta$NSC,method = 'pearson')
# corr_NSC_5 = cor(umaps$X5,meta$NSC,method = 'pearson')
# corr_NSC_6 = cor(umaps$X6,meta$NSC,method = 'pearson')
# 
# corr_RSC_1 = cor(umaps$X1,meta$RSC,method = 'pearson')
# corr_RSC_2 = cor(umaps$X2,meta$RSC,method = 'pearson')
# corr_RSC_3 = cor(umaps$X3,meta$RSC,method = 'pearson')
# corr_RSC_4 = cor(umaps$X4,meta$RSC,method = 'pearson')
# corr_RSC_5 = cor(umaps$X5,meta$RSC,method = 'pearson')
# corr_RSC_6 = cor(umaps$X6,meta$RSC,method = 'pearson')
# 
# corr_FRiP_1 = cor(umaps$X1,meta$FRiP,method = 'pearson')
# corr_FRiP_2 = cor(umaps$X2,meta$FRiP,method = 'pearson')
# corr_FRiP_3 = cor(umaps$X3,meta$FRiP,method = 'pearson')
# corr_FRiP_4 = cor(umaps$X4,meta$FRiP,method = 'pearson')
# corr_FRiP_5 = cor(umaps$X5,meta$FRiP,method = 'pearson')
# corr_FRiP_6 = cor(umaps$X6,meta$FRiP,method = 'pearson')
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
# group_list2 = meta$unique_reads
# group_list3 = meta$Peaks
# group_list4 = meta$NSC
# group_list5 = meta$RSC
# group_list6 = meta$FRiP
# group_list7 = meta$sequencing_depth
# group_list8 = meta$batch
# 
# colData=data.frame(row.names = colnames(group_D1A1),
#                    condition=group_list1,
#                    reads=group_list2,
#                    peaks=group_list3,
#                    NSC=group_list4,
#                    RSC = group_list5,
#                    FRiP = group_list6,
#                    sequencing_depth = group_list7,
#                    batch = group_list8)
# 
# colData$condition = factor(colData$condition,levels = c('D','A'))
# colData$reads =as.numeric(colData$reads)
# colData$peaks=as.numeric(colData$peaks)
# colData$NSC=as.numeric(colData$NSC)
# colData$RSC=as.numeric(colData$RSC)
# colData$FRiP=as.numeric(colData$FRiP)
# colData$sequencing_depth=as.numeric(colData$sequencing_depth)
# colData$batch=factor(colData$batch,)
# 
# dds=DESeqDataSetFromMatrix(countData = group_D1A1,
#                            colData = colData,
#                            design = ~ condition )
# #tmp = model.matrix(~group_list + gender + Age + PMI, colData)
# 
# dds=DESeq(dds)
# 
# res_glia=results(dds,
#                  contrast = c("condition","D","A"))
# resOrdered_glia=res_glia[order(res_glia$padj),]
# DEG_glia=as.data.frame(resOrdered_glia)
# DEG_glia=na.omit(DEG_glia)
# nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.6),]
# # nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.5),]
# choose_gene_neuron=rownames(nrDEG_glia)
# choose_matrix=group_D1A1[choose_gene_neuron,]
# DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
# up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0.5),]
# down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < -0.5),]
# diff_peak = peak_info[rownames(DEG_matrix),]
# new_DEG_matrix = bind_cols(diff_peak,DEG_matrix)
# setwd("~/Documents/Virginia Tech/Project_3_MIA mouse/differential analysis/results/new_promoters/ABCD/differential peaks/D1A1")
# all_glia = group_D1A1[rownames(DEG_glia),]
# all_matrix = data.frame(DEG_glia,all_glia)
# write.csv(all_matrix,'All_D1A1.csv')
# write.csv(new_DEG_matrix,'Differential_D1A1.csv')
# 
# ###############################D_2 vs. A_2###############################################################
# meta = bind_rows(factors[66:71,1:8],factors[12:17,1:8])
# meta$condition_num = c(rep(1,6),rep(2,6))
# group_D2A2 = data.frame(def[,66:71],def[,12:17])
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
# cor_con_peaks = cor(meta$condition_num,meta$Peaks,method = 'pearson')
# cor_con_NSC = cor(meta$condition_num,meta$NSC,method = 'pearson')
# cor_con_RSC = cor(meta$condition_num,meta$RSC,method = 'pearson')
# cor_con_FRiP = cor(meta$condition_num,meta$FRiP,method = 'pearson')
# cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
# cor_con_uni_reads = cor(meta$condition_num,meta$unique_reads,method = 'pearson')
# # FRiP NSC 
# # FRiP NSC -UMAP
# corr_uni_reads_1 = cor(umaps$X1,meta$unique_reads,method = 'pearson')
# corr_uni_reads_2 = cor(umaps$X2,meta$unique_reads,method = 'pearson')
# corr_uni_reads_3 = cor(umaps$X3,meta$unique_reads,method = 'pearson')
# corr_uni_reads_4 = cor(umaps$X4,meta$unique_reads,method = 'pearson')
# corr_uni_reads_5 = cor(umaps$X5,meta$unique_reads,method = 'pearson')
# corr_uni_reads_6 = cor(umaps$X6,meta$unique_reads,method = 'pearson')
# 
# corr_peaks_1 = cor(umaps$X1,meta$Peaks,method = 'pearson')
# corr_peaks_2 = cor(umaps$X2,meta$Peaks,method = 'pearson')
# corr_peaks_3 = cor(umaps$X3,meta$Peaks,method = 'pearson')
# corr_peaks_4 = cor(umaps$X4,meta$Peaks,method = 'pearson')
# corr_peaks_5 = cor(umaps$X5,meta$Peaks,method = 'pearson')
# corr_peaks_6 = cor(umaps$X6,meta$Peaks,method = 'pearson')
# 
# corr_NSC_1 = cor(umaps$X1,meta$NSC,method = 'pearson')
# corr_NSC_2 = cor(umaps$X2,meta$NSC,method = 'pearson')
# corr_NSC_3 = cor(umaps$X3,meta$NSC,method = 'pearson')
# corr_NSC_4 = cor(umaps$X4,meta$NSC,method = 'pearson')
# corr_NSC_5 = cor(umaps$X5,meta$NSC,method = 'pearson')
# corr_NSC_6 = cor(umaps$X6,meta$NSC,method = 'pearson')
# 
# corr_RSC_1 = cor(umaps$X1,meta$RSC,method = 'pearson')
# corr_RSC_2 = cor(umaps$X2,meta$RSC,method = 'pearson')
# corr_RSC_3 = cor(umaps$X3,meta$RSC,method = 'pearson')
# corr_RSC_4 = cor(umaps$X4,meta$RSC,method = 'pearson')
# corr_RSC_5 = cor(umaps$X5,meta$RSC,method = 'pearson')
# corr_RSC_6 = cor(umaps$X6,meta$RSC,method = 'pearson')
# 
# corr_FRiP_1 = cor(umaps$X1,meta$FRiP,method = 'pearson')
# corr_FRiP_2 = cor(umaps$X2,meta$FRiP,method = 'pearson')
# corr_FRiP_3 = cor(umaps$X3,meta$FRiP,method = 'pearson')
# corr_FRiP_4 = cor(umaps$X4,meta$FRiP,method = 'pearson')
# corr_FRiP_5 = cor(umaps$X5,meta$FRiP,method = 'pearson')
# corr_FRiP_6 = cor(umaps$X6,meta$FRiP,method = 'pearson')
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
# group_list2 = meta$unique_reads
# group_list3 = meta$Peaks
# group_list4 = meta$NSC
# group_list5 = meta$RSC
# group_list6 = meta$FRiP
# group_list7 = meta$sequencing_depth
# group_list8 = meta$batch
# 
# colData=data.frame(row.names = colnames(group_D2A2),
#                    condition=group_list1,
#                    reads=group_list2,
#                    peaks=group_list3,
#                    NSC=group_list4,
#                    RSC = group_list5,
#                    FRiP = group_list6,
#                    sequencing_depth = group_list7,
#                    batch = group_list8)
# 
# colData$condition = factor(colData$condition,levels = c('D','A'))
# colData$reads =as.numeric(colData$reads)
# colData$peaks=as.numeric(colData$peaks)
# colData$NSC=as.numeric(colData$NSC)
# colData$RSC=as.numeric(colData$RSC)
# colData$batch=as.factor(colData$batch)
# colData$FRiP=as.numeric(colData$FRiP)
# colData$sequencing_depth=as.numeric(colData$sequencing_depth)
# 
# dds=DESeqDataSetFromMatrix(countData = group_D2A2,
#                            colData = colData,
#                            design = ~ condition + peaks)
# #tmp = model.matrix(~group_list + gender + Age + PMI, colData)
# 
# dds=DESeq(dds)
# 
# res_glia=results(dds,
#                  contrast = c("condition","D","A"))
# resOrdered_glia=res_glia[order(res_glia$padj),]
# DEG_glia=as.data.frame(resOrdered_glia)
# DEG_glia=na.omit(DEG_glia)
# nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.6),]
# # nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.5),]
# choose_gene_neuron=rownames(nrDEG_glia)
# choose_matrix=group_D2A2[choose_gene_neuron,]
# DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
# up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0 & nrDEG_glia$log2FoldChange > 0.5),]
# down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0 & nrDEG_glia$log2FoldChange < -0.5),]
# diff_peak = peak_info[rownames(DEG_matrix),]
# new_DEG_matrix = bind_cols(diff_peak,DEG_matrix)
# setwd("~/Documents/Virginia Tech/Project_3_MIA mouse/differential analysis/results/new_promoters/ABCD/differential peaks/D2A2")
# all_glia = group_D2A2[rownames(DEG_glia),]
# all_matrix = data.frame(DEG_glia,all_glia)
# write.csv(all_matrix,'All_D2A2.csv')
# write.csv(new_DEG_matrix,'Differential_D2A2.csv')
# ###############################D_3 vs. A_3###############################################################
# meta = bind_rows(factors[72:82,1:8],factors[18:29,1:8])
# meta$condition_num = c(rep(1,11),rep(2,12))
# group_D3A3 = data.frame(def[,72:82],def[,18:29])
# miss_1 <- c()
# for (i in 1:nrow(group_D3A3)) {
#   if(length(which(group_D3A3[i,1:11]<10)) > 0.5*ncol(group_D3A3[,1:11])) miss_1 <- append(miss_1,i)
# }
# miss_2 <- c()
# for (i in 1:nrow(group_D3A3)) {
#   if(length(which(group_D3A3[i,12:23]<10)) > 0.5*ncol(group_D3A3[,12:23])) miss_2 <- append(miss_2,i)
# }
# miss = union(miss_1,miss_2)
# group_D3A3 = group_D3A3[-miss,]
# t_group_D3A3 = t(group_D3A3)
# umap_results <- umap(t_group_D3A3,n_neighbors=11,min_dist=0.1,n_components = 6)
# umaps = data.frame(umap_results[["layout"]])
# umaps$condition = as.factor(meta$condition)
# ggplot(umaps,mapping = aes(x=X1,y=X2,col=condition)) +
#   geom_point(aes(colour = condition))
# 
# cor_con_batch = cor(meta$condition_num,meta$batch,method = 'pearson')
# cor_con_peaks = cor(meta$condition_num,meta$Peaks,method = 'pearson')
# cor_con_NSC = cor(meta$condition_num,meta$NSC,method = 'pearson')
# cor_con_RSC = cor(meta$condition_num,meta$RSC,method = 'pearson')
# cor_con_FRiP = cor(meta$condition_num,meta$FRiP,method = 'pearson')
# cor_con_sequencing_depth = cor(meta$condition_num,meta$sequencing_depth,method = 'pearson')
# cor_con_uni_reads = cor(meta$condition_num,meta$unique_reads,method = 'pearson')
# # 
# # RSC reads -UMAP
# corr_uni_reads_1 = cor(umaps$X1,meta$unique_reads,method = 'pearson')
# corr_uni_reads_2 = cor(umaps$X2,meta$unique_reads,method = 'pearson')
# corr_uni_reads_3 = cor(umaps$X3,meta$unique_reads,method = 'pearson')
# corr_uni_reads_4 = cor(umaps$X4,meta$unique_reads,method = 'pearson')
# corr_uni_reads_5 = cor(umaps$X5,meta$unique_reads,method = 'pearson')
# corr_uni_reads_6 = cor(umaps$X6,meta$unique_reads,method = 'pearson')
# 
# corr_peaks_1 = cor(umaps$X1,meta$Peaks,method = 'pearson')
# corr_peaks_2 = cor(umaps$X2,meta$Peaks,method = 'pearson')
# corr_peaks_3 = cor(umaps$X3,meta$Peaks,method = 'pearson')
# corr_peaks_4 = cor(umaps$X4,meta$Peaks,method = 'pearson')
# corr_peaks_5 = cor(umaps$X5,meta$Peaks,method = 'pearson')
# corr_peaks_6 = cor(umaps$X6,meta$Peaks,method = 'pearson')
# 
# corr_NSC_1 = cor(umaps$X1,meta$NSC,method = 'pearson')
# corr_NSC_2 = cor(umaps$X2,meta$NSC,method = 'pearson')
# corr_NSC_3 = cor(umaps$X3,meta$NSC,method = 'pearson')
# corr_NSC_4 = cor(umaps$X4,meta$NSC,method = 'pearson')
# corr_NSC_5 = cor(umaps$X5,meta$NSC,method = 'pearson')
# corr_NSC_6 = cor(umaps$X6,meta$NSC,method = 'pearson')
# 
# corr_RSC_1 = cor(umaps$X1,meta$RSC,method = 'pearson')
# corr_RSC_2 = cor(umaps$X2,meta$RSC,method = 'pearson')
# corr_RSC_3 = cor(umaps$X3,meta$RSC,method = 'pearson')
# corr_RSC_4 = cor(umaps$X4,meta$RSC,method = 'pearson')
# corr_RSC_5 = cor(umaps$X5,meta$RSC,method = 'pearson')
# corr_RSC_6 = cor(umaps$X6,meta$RSC,method = 'pearson')
# 
# corr_FRiP_1 = cor(umaps$X1,meta$FRiP,method = 'pearson')
# corr_FRiP_2 = cor(umaps$X2,meta$FRiP,method = 'pearson')
# corr_FRiP_3 = cor(umaps$X3,meta$FRiP,method = 'pearson')
# corr_FRiP_4 = cor(umaps$X4,meta$FRiP,method = 'pearson')
# corr_FRiP_5 = cor(umaps$X5,meta$FRiP,method = 'pearson')
# corr_FRiP_6 = cor(umaps$X6,meta$FRiP,method = 'pearson')
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
# group_list2 = meta$unique_reads
# group_list3 = meta$Peaks
# group_list4 = meta$NSC
# group_list5 = meta$RSC
# group_list6 = meta$FRiP
# group_list7 = meta$sequencing_depth
# group_list8 = meta$batch
# 
# colData=data.frame(row.names = colnames(group_D3A3),
#                    condition=group_list1,
#                    reads=group_list2,
#                    peaks=group_list3,
#                    NSC=group_list4,
#                    RSC = group_list5,
#                    FRiP = group_list6,
#                    sequencing_depth = group_list7,
#                    batch = group_list8)
# 
# colData$condition = factor(colData$condition,levels = c('D','A'))
# colData$reads =as.numeric(colData$reads)
# colData$peaks=as.numeric(colData$peaks)
# colData$NSC=as.numeric(colData$NSC)
# colData$RSC=as.numeric(colData$RSC)
# colData$FRiP=as.numeric(colData$FRiP)
# colData$batch=as.factor(colData$batch)
# colData$sequencing_depth=as.numeric(colData$sequencing_depth)
# 
# dds=DESeqDataSetFromMatrix(countData = group_D3A3,
#                            colData = colData,
#                            design = ~ condition + batch + FRiP + peaks + RSC + reads)
# #tmp = model.matrix(~group_list + gender + Age + PMI, colData)
# 
# dds=DESeq(dds)
# 
# res_glia=results(dds,
#                  contrast = c("condition","D","A"))
# resOrdered_glia=res_glia[order(res_glia$padj),]
# DEG_glia=as.data.frame(resOrdered_glia)
# DEG_glia=na.omit(DEG_glia)
# nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.6),]
# # nrDEG_glia=DEG_glia[which(DEG_glia$padj<0.05 & abs(DEG_glia$log2FoldChange) > 0.5),]
# choose_gene_neuron=rownames(nrDEG_glia)
# choose_matrix=group_D3A3[choose_gene_neuron,]
# DEG_matrix = data.frame(nrDEG_glia,choose_matrix)
# up = nrDEG_glia[which(nrDEG_glia$log2FoldChange > 0 & nrDEG_glia$log2FoldChange > 0.5),]
# down = nrDEG_glia[which(nrDEG_glia$log2FoldChange < 0 & nrDEG_glia$log2FoldChange < -0.5),]
# diff_peak = peak_info[rownames(DEG_matrix),]
# new_DEG_matrix = bind_cols(diff_peak,DEG_matrix)
# setwd("~/Documents/Virginia Tech/Project_3_MIA mouse/differential analysis/results/new_promoters/ABCD/differential peaks/D3A3")
# all_glia = group_D3A3[rownames(DEG_glia),]
# all_matrix = data.frame(DEG_glia,all_glia)
# write.csv(all_matrix,'All_D3A3.csv')
# write.csv(new_DEG_matrix,'Differential_D3A3.csv')

