library(ChIPseeker)
#BiocManager::install("TxDb.Mmusculus.UCSC.mm9.knownGene")
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

#library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(stringr)
library(reshape2)
library(ggplot2)
library(ggpubr)

rm(list = ls())
options(stringsAsFactors = F)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
##############H3K4me3_ABCD###############
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/count_files/")
bedPeaksFile = "H3K4me3_peak_abcd.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
cat(paste0('There are ', length(peak), 'peaks for this data'))
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
anno = grepl('Promoter', peakAnno_df$annotation,fixed = TRUE)
new_peakAnno_df = peakAnno_df[anno,]

H3K4me3_peak = read.table('H3K4me3_peak_abcd.bed')
H3K4me3 = read.csv('H3K4me3_count_abcd.csv')
colnames(H3K4me3_peak) = colnames(H3K4me3)[1:4]
promoter_peak = H3K4me3_peak[which(H3K4me3_peak$peak_num %in% new_peakAnno_df$V4),]
promoter_count = H3K4me3[which(H3K4me3$peak_num %in% new_peakAnno_df$V4),]
write.csv(promoter_count,'promoter_abcd_count.csv',row.names = FALSE)
write.table(promoter_peak,'promoter_abcd.bed',
            row.names = FALSE, sep="\t", quote = FALSE)
# #####################raw_count########################################
# bedPeaksFile = "H3K4me3_peak_abcd_raw.bed"
# peak <- readPeakFile(bedPeaksFile)
# keepChr=!grepl('_',seqlevels(peak))
# seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
# cat(paste0('There are ', length(peak), 'peaks for this data'))
# peakAnno <- annotatePeak(peak,
#                          tssRegion = c(-2000, 2000),
#                          TxDb = txdb,
#                          annoDb = "org.Mm.eg.db")
# peakAnno_df = as.data.frame(peakAnno)
# anno = grepl('Promoter', peakAnno_df$annotation,fixed = TRUE)
# new_peakAnno_df = peakAnno_df[anno,]
# 
# H3K4me3_peak = read.table('H3K4me3_peak_abcd_raw.bed')
# H3K4me3 = read.csv('H3K4me3_count_abcd_raw.csv')
# colnames(H3K4me3_peak) = colnames(H3K4me3)[1:4]
# promoter_peak = H3K4me3_peak[which(H3K4me3_peak$peak_num %in% new_peakAnno_df$V4),]
# promoter_count = H3K4me3[which(H3K4me3$peak_num %in% new_peakAnno_df$V4),]
# write.csv(promoter_count,'promoter_abcd_count_raw.csv',row.names = FALSE)
# write.table(promoter_peak,'promoter_abcd_raw.bed',
#             row.names = FALSE, sep="\t", quote = FALSE)
######################H3K27ac_ABCD####################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/count_files/")
bedPeaksFile = "H3K27ac_peak_abcd.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
cat(paste0('There are ', length(peak), 'peaks for this data'))
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
anno = grepl('Promoter', peakAnno_df$annotation,fixed = TRUE)
new_peakAnno_df = peakAnno_df[!anno,]

H3K27ac_peak = read.table('H3K27ac_peak_abcd.bed')
H3K27ac = read.csv('H3K27ac_count_abcd.csv')
colnames(H3K27ac_peak) = colnames(H3K27ac)[1:4]
enhancer_peak = H3K27ac_peak[which(H3K27ac_peak$peak_num %in% new_peakAnno_df$V4),]
enhancer_count = H3K27ac[which(H3K27ac$peak_num %in% new_peakAnno_df$V4),]
write.csv(enhancer_count,'enhancer_count_abcd.csv',row.names = FALSE)
write.table(enhancer_peak,'enhancer_abcd.bed',
            row.names = FALSE, sep="\t", quote = FALSE)

# #######################raw_count######################################
# bedPeaksFile = "H3K27ac_peak_abcd_raw.bed"
# peak <- readPeakFile(bedPeaksFile)
# keepChr=!grepl('_',seqlevels(peak))
# seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
# cat(paste0('There are ', length(peak), 'peaks for this data'))
# peakAnno <- annotatePeak(peak,
#                          tssRegion = c(-2000, 2000),
#                          TxDb = txdb,
#                          annoDb = "org.Mm.eg.db")
# peakAnno_df = as.data.frame(peakAnno)
# anno = grepl('Promoter', peakAnno_df$annotation,fixed = TRUE)
# new_peakAnno_df = peakAnno_df[!anno,]
# 
# H3K27ac_peak = read.table('H3K27ac_peak_abcd_raw.bed')
# H3K27ac = read.csv('H3K27ac_count_abcd_raw.csv')
# colnames(H3K27ac_peak) = colnames(H3K27ac)[1:4]
# enhancer_peak = H3K27ac_peak[which(H3K27ac_peak$peak_num %in% new_peakAnno_df$V4),]
# enhancer_count = H3K27ac[which(H3K27ac$peak_num %in% new_peakAnno_df$V4),]
# write.csv(enhancer_count,'enhancer_count_abcd_raw.csv',row.names = FALSE)
# write.table(enhancer_peak,'enhancer_abcd_raw.bed',
#             row.names = FALSE, sep="\t", quote = FALSE)
# 
##############H3K4me3_EFG###############
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/count_files/")
bedPeaksFile = "H3K4me3_peak_efg.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
cat(paste0('There are ', length(peak), 'peaks for this data'))
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
anno = grepl('Promoter', peakAnno_df$annotation,fixed = TRUE)
new_peakAnno_df = peakAnno_df[anno,]

H3K4me3_peak = read.table('H3K4me3_peak_efg.bed')
H3K4me3 = read.csv('H3K4me3_count_efg.csv')
colnames(H3K4me3_peak) = colnames(H3K4me3)[1:4]
promoter_peak = H3K4me3_peak[which(H3K4me3_peak$peak_num %in% new_peakAnno_df$V4),]
promoter_count = H3K4me3[which(H3K4me3$peak_num %in% new_peakAnno_df$V4),]
write.csv(promoter_count,'promoter_efg_count.csv',row.names = FALSE)
write.table(promoter_peak,'promoter_efg.bed',
            row.names = FALSE, sep="\t", quote = FALSE)


# ##############################raw_count####################################
# bedPeaksFile = "H3K4me3_peak_efg_raw.bed"
# peak <- readPeakFile(bedPeaksFile)
# keepChr=!grepl('_',seqlevels(peak))
# seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
# cat(paste0('There are ', length(peak), 'peaks for this data'))
# peakAnno <- annotatePeak(peak,
#                          tssRegion = c(-2000, 2000),
#                          TxDb = txdb,
#                          annoDb = "org.Mm.eg.db")
# peakAnno_df = as.data.frame(peakAnno)
# anno = grepl('Promoter', peakAnno_df$annotation,fixed = TRUE)
# new_peakAnno_df = peakAnno_df[anno,]
# 
# H3K4me3_peak = read.table('H3K4me3_peak_efg_raw.bed')
# H3K4me3 = read.csv('H3K4me3_count_efg_raw.csv')
# colnames(H3K4me3_peak) = colnames(H3K4me3)[1:4]
# promoter_peak = H3K4me3_peak[which(H3K4me3_peak$peak_num %in% new_peakAnno_df$V4),]
# promoter_count = H3K4me3[which(H3K4me3$peak_num %in% new_peakAnno_df$V4),]
# write.csv(promoter_count,'promoter_efg_count_raw.csv',row.names = FALSE)
# write.table(promoter_peak,'promoter_efg_raw.bed',
#             row.names = FALSE, sep="\t", quote = FALSE)
######################H3K27ac_EFG####################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/count_files/")
bedPeaksFile = "H3K27ac_peak_efg.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
cat(paste0('There are ', length(peak), 'peaks for this data'))
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
anno = grepl('Promoter', peakAnno_df$annotation,fixed = TRUE)
new_peakAnno_df = peakAnno_df[!anno,]

H3K27ac_peak = read.table('H3K27ac_peak_efg.bed')
H3K27ac = read.csv('H3K27ac_count_efg.csv')
colnames(H3K27ac_peak) = colnames(H3K27ac)[1:4]
enhancer_peak = H3K27ac_peak[which(H3K27ac_peak$peak_num %in% new_peakAnno_df$V4),]
enhancer_count = H3K27ac[which(H3K27ac$peak_num %in% new_peakAnno_df$V4),]
write.csv(enhancer_count,'enhancer_count_efg.csv',row.names = FALSE)
write.table(enhancer_peak,'enhancer_efg.bed',
            row.names = FALSE, sep="\t", quote = FALSE)

# ######################raw_count###################################
# 
# bedPeaksFile = "H3K27ac_peak_efg_raw.bed"
# peak <- readPeakFile(bedPeaksFile)
# keepChr=!grepl('_',seqlevels(peak))
# seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
# cat(paste0('There are ', length(peak), 'peaks for this data'))
# peakAnno <- annotatePeak(peak,
#                          tssRegion = c(-2000, 2000),
#                          TxDb = txdb,
#                          annoDb = "org.Mm.eg.db")
# peakAnno_df = as.data.frame(peakAnno)
# anno = grepl('Promoter', peakAnno_df$annotation,fixed = TRUE)
# new_peakAnno_df = peakAnno_df[!anno,]
# 
# H3K27ac_peak = read.table('H3K27ac_peak_efg_raw.bed')
# H3K27ac = read.csv('H3K27ac_count_efg_raw.csv')
# colnames(H3K27ac_peak) = colnames(H3K27ac)[1:4]
# enhancer_peak = H3K27ac_peak[which(H3K27ac_peak$peak_num %in% new_peakAnno_df$V4),]
# enhancer_count = H3K27ac[which(H3K27ac$peak_num %in% new_peakAnno_df$V4),]
# write.csv(enhancer_count,'enhancer_count_efg_raw.csv',row.names = FALSE)
# write.table(enhancer_peak,'enhancer_efg_raw.bed',
#             row.names = FALSE, sep="\t", quote = FALSE)