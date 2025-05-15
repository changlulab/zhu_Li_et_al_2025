rm(list = ls())
options(stringsAsFactors = F)
library(DiffBind)
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/diffbind/inputs_sheets")
########################################################################################################
H3K4me3 <- dba(sampleSheet = "H3K4me3_abcd.csv")
H3K4me3_blacklist_remove = dba.blacklist(H3K4me3,blacklist = DBA_BLACKLIST_MM10, greylist = FALSE)
H3K4me3_consensus <- dba.peakset(H3K4me3_blacklist_remove,
                                 consensus = DBA_TREATMENT, 
                                 minOverlap = 2)
H3K4me3_consensus
H3K4me3_consensus_set <- dba(H3K4me3_consensus,
                             mask = H3K4me3_consensus$masks$Consensus, 
                             minOverlap = 1)
H3K4me3_consensus_set
consensus_peaks <- dba.peakset(H3K4me3_consensus_set, bRetrieve = TRUE)
H3K4me3_count <- dba.count(H3K4me3_blacklist_remove, 
                           summits = FALSE,
                           peaks = consensus_peaks,filter=1,
                           bScaleControl = TRUE,
                           minCount=1,
                           score=DBA_SCORE_TMM_MINUS_FULL)

setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result")
save(H3K4me3,H3K4me3_blacklist_remove,H3K4me3_consensus,
     H3K4me3_consensus_set,consensus_peaks,H3K4me3_count,
     file = 'H3K4me3_abcd.RData')

H3K4me3_count_def = data.frame(H3K4me3_count[["binding"]])
H3K4me3_count[["chrmap"]]
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '1')] = 'chr1'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '2')] = 'chr10'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '3')] = 'chr11'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '4')] = 'chr12'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '5')] = 'chr13'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '6')] = 'chr14'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '7')] = 'chr15'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '8')] = 'chr16'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '9')] = 'chr17'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '10')] = 'chr18'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '11')] = 'chr19'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '12')] = 'chr2'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '13')] = 'chr3'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '14')] = 'chr4'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '15')] = 'chr5'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '16')] = 'chr6'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '17')] = 'chr7'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '18')] = 'chr8'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '19')] = 'chr9'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '20')] = 'chrX'
H3K4me3_count_def$CHR[which(H3K4me3_count_def$CHR == '21')] = 'chrY'
chr_list = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
             "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
             "chrX","chrY")

miss <- c()
for (i in 1:nrow(H3K4me3_count_def)) {
  if(!(H3K4me3_count_def$CHR[i] %in% chr_list)) miss <- append(miss,i)
}
# new_H3K4me3_count_def = H3K4me3_count_def[-miss,]
H3K4me3_count_def[,4:ncol(H3K4me3_count_def)] = round(H3K4me3_count_def[,4:ncol(H3K4me3_count_def)])
new_H3K4me3_count_def = data.frame(H3K4me3_count_def[,1:3],
                                       rownames(H3K4me3_count_def),
                                       H3K4me3_count_def[,4:ncol(H3K4me3_count_def)])
colnames(new_H3K4me3_count_def)[1:4] = c('chr','start','end','peak_num')
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/count_files")
write.csv(new_H3K4me3_count_def,'H3K4me3_count_abcd.csv',row.names = FALSE)
write.table(new_H3K4me3_count_def[,1:4],'H3K4me3_peak_abcd.bed',sep = "\t",
            quote = F,
            col.names = F,
            row.names = F)

# ##################raw_count###############################################
# H3K4me3_count_raw <- dba.count(H3K4me3_blacklist_remove, 
#                                summits = FALSE,
#                                peaks = consensus_peaks,filter=1,
#                                bScaleControl = TRUE,
#                                minCount=1,
#                                score=DBA_SCORE_READS)
# 
# H3K4me3_count_raw_def = data.frame(H3K4me3_count_raw[["binding"]])
# H3K4me3_count_raw[["chrmap"]]
# H3K4me3_count_raw_def$CHR[which(H3K4me3_count_raw_def$CHR == '1')] = 'chr1'
# H3K4me3_count_raw_def$CHR[which(H3K4me3_count_raw_def$CHR == '2')] = 'chr10'
# H3K4me3_count_raw_def$CHR[which(H3K4me3_count_raw_def$CHR == '3')] = 'chr11'
# H3K4me3_count_raw_def$CHR[which(H3K4me3_count_raw_def$CHR == '4')] = 'chr12'
# H3K4me3_count_raw_def$CHR[which(H3K4me3_count_raw_def$CHR == '5')] = 'chr13'
# H3K4me3_count_raw_def$CHR[which(H3K4me3_count_raw_def$CHR == '6')] = 'chr14'
# H3K4me3_count_raw_def$CHR[which(H3K4me3_count_raw_def$CHR == '7')] = 'chr15'
# H3K4me3_count_raw_def$CHR[which(H3K4me3_count_raw_def$CHR == '8')] = 'chr16'
# H3K4me3_count_raw_def$CHR[which(H3K4me3_count_raw_def$CHR == '9')] = 'chr17'
# H3K4me3_count_raw_def$CHR[which(H3K4me3_count_raw_def$CHR == '10')] = 'chr18'
# H3K4me3_count_raw_def$CHR[which(H3K4me3_count_raw_def$CHR == '11')] = 'chr19'
# H3K4me3_count_raw_def$CHR[which(H3K4me3_count_raw_def$CHR == '12')] = 'chr2'
# H3K4me3_count_raw_def$CHR[which(H3K4me3_count_raw_def$CHR == '13')] = 'chr3'
# H3K4me3_count_raw_def$CHR[which(H3K4me3_count_raw_def$CHR == '14')] = 'chr4'
# H3K4me3_count_raw_def$CHR[which(H3K4me3_count_raw_def$CHR == '15')] = 'chr5'
# H3K4me3_count_raw_def$CHR[which(H3K4me3_count_raw_def$CHR == '16')] = 'chr6'
# H3K4me3_count_raw_def$CHR[which(H3K4me3_count_raw_def$CHR == '17')] = 'chr7'
# H3K4me3_count_raw_def$CHR[which(H3K4me3_count_raw_def$CHR == '18')] = 'chr8'
# H3K4me3_count_raw_def$CHR[which(H3K4me3_count_raw_def$CHR == '19')] = 'chr9'
# H3K4me3_count_raw_def$CHR[which(H3K4me3_count_raw_def$CHR == '20')] = 'chrX'
# H3K4me3_count_raw_def$CHR[which(H3K4me3_count_raw_def$CHR == '21')] = 'chrY'
# chr_list = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
#              "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
#              "chrX","chrY")
# miss <- c()
# for (i in 1:nrow(H3K4me3_count_raw_def)) {
#   if(!(H3K4me3_count_raw_def$CHR[i] %in% chr_list)) miss <- append(miss,i)
# }
# 
# H3K4me3_count_raw_def[,4:ncol(H3K4me3_count_raw_def)] = round(H3K4me3_count_raw_def[,4:ncol(H3K4me3_count_raw_def)])
# new_H3K4me3_count_raw_def = data.frame(H3K4me3_count_raw_def[,1:3],
#                                    rownames(H3K4me3_count_raw_def),
#                                    H3K4me3_count_raw_def[,4:ncol(H3K4me3_count_raw_def)])
# colnames(new_H3K4me3_count_raw_def)[1:4] = c('chr','start','end','peak_num')
# setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/count_files")
# write.csv(new_H3K4me3_count_raw_def,'H3K4me3_count_abcd_raw.csv',row.names = FALSE)
# write.table(new_H3K4me3_count_raw_def[,1:4],'H3K4me3_peak_abcd_raw.bed',sep = "\t",
#             quote = F,
#             col.names = F,
#             row.names = F)

#########################################################################################################
rm(list = ls())
library(DiffBind)
options(stringsAsFactors = F)
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/diffbind/inputs_sheets")
H3K27ac <- dba(sampleSheet = "H3K27ac_abcd.csv")
H3K27ac_blacklist_remove = dba.blacklist(H3K27ac,blacklist = DBA_BLACKLIST_MM10, greylist = TRUE)
H3K27ac_consensus <- dba.peakset(H3K27ac_blacklist_remove,
                                 consensus = DBA_TREATMENT, 
                                 minOverlap = 2)
H3K27ac_consensus_set <- dba(H3K27ac_consensus,
                             mask = H3K27ac_consensus$masks$Consensus, 
                             minOverlap = 1)
consensus_peaks <- dba.peakset(H3K27ac_consensus_set, bRetrieve = TRUE)
H3K27ac_count <- dba.count(H3K27ac_blacklist_remove, 
                           summits = FALSE,
                           peaks = consensus_peaks,filter=1,
                           bScaleControl = TRUE,
                           minCount=1,
                           score=DBA_SCORE_TMM_MINUS_FULL)
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result")
save(H3K27ac,H3K27ac_blacklist_remove,H3K27ac_consensus,
     H3K27ac_consensus_set,consensus_peaks,H3K27ac_count,
     file = 'H3K27ac_abcd.RData')

H3K27ac_count_def = data.frame(H3K27ac_count[["binding"]])
H3K27ac_count[["chrmap"]]
H3K27ac_count_def$CHR[which(H3K27ac_count_def$CHR == '1')] = 'chr1'
H3K27ac_count_def$CHR[which(H3K27ac_count_def$CHR == '2')] = 'chr10'
H3K27ac_count_def$CHR[which(H3K27ac_count_def$CHR == '3')] = 'chr11'
H3K27ac_count_def$CHR[which(H3K27ac_count_def$CHR == '4')] = 'chr12'
H3K27ac_count_def$CHR[which(H3K27ac_count_def$CHR == '5')] = 'chr13'
H3K27ac_count_def$CHR[which(H3K27ac_count_def$CHR == '6')] = 'chr14'
H3K27ac_count_def$CHR[which(H3K27ac_count_def$CHR == '7')] = 'chr15'
H3K27ac_count_def$CHR[which(H3K27ac_count_def$CHR == '8')] = 'chr16'
H3K27ac_count_def$CHR[which(H3K27ac_count_def$CHR == '9')] = 'chr17'
H3K27ac_count_def$CHR[which(H3K27ac_count_def$CHR == '10')] = 'chr18'
H3K27ac_count_def$CHR[which(H3K27ac_count_def$CHR == '11')] = 'chr19'
H3K27ac_count_def$CHR[which(H3K27ac_count_def$CHR == '12')] = 'chr2'
H3K27ac_count_def$CHR[which(H3K27ac_count_def$CHR == '13')] = 'chr3'
H3K27ac_count_def$CHR[which(H3K27ac_count_def$CHR == '14')] = 'chr4'
H3K27ac_count_def$CHR[which(H3K27ac_count_def$CHR == '15')] = 'chr5'
H3K27ac_count_def$CHR[which(H3K27ac_count_def$CHR == '16')] = 'chr6'
H3K27ac_count_def$CHR[which(H3K27ac_count_def$CHR == '17')] = 'chr7'
H3K27ac_count_def$CHR[which(H3K27ac_count_def$CHR == '18')] = 'chr8'
H3K27ac_count_def$CHR[which(H3K27ac_count_def$CHR == '19')] = 'chr9'
H3K27ac_count_def$CHR[which(H3K27ac_count_def$CHR == '20')] = 'chrX'
H3K27ac_count_def$CHR[which(H3K27ac_count_def$CHR == '21')] = 'chrY'
chr_list = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
             "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
             "chrX","chrY")
miss <- c()
for (i in 1:nrow(H3K27ac_count_def)) {
  if(!(H3K27ac_count_def$CHR[i] %in% chr_list)) miss <- append(miss,i)
}
# new_H3K27ac_count_def = H3K27ac_count_def[-miss,]
H3K27ac_count_def[,4:ncol(H3K27ac_count_def)] = round(H3K27ac_count_def[,4:ncol(H3K27ac_count_def)])
new_H3K27ac_count_def = data.frame(H3K27ac_count_def[,1:3],
                                   rownames(H3K27ac_count_def),
                                   H3K27ac_count_def[,4:ncol(H3K27ac_count_def)])
colnames(new_H3K27ac_count_def)[1:4] = c('chr','start','end','peak_num')
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/count_files")
write.csv(new_H3K27ac_count_def,'H3K27ac_count_abcd.csv',row.names = FALSE)
write.table(new_H3K27ac_count_def[,1:4],'H3K27ac_peak_abcd.bed',sep = "\t",
            quote = F,
            col.names = F,
            row.names = F)

# ##################raw_count###############################################
# H3K27ac_count_raw <- dba.count(H3K27ac_blacklist_remove, 
#                            summits = FALSE,
#                            peaks = consensus_peaks,filter=1,
#                            bScaleControl = TRUE,
#                            minCount=1,
#                            score=DBA_SCORE_READS)
# 
# 
# H3K27ac_count_raw_def = data.frame(H3K27ac_count_raw[["binding"]])
# H3K27ac_count_raw[["chrmap"]]
# H3K27ac_count_raw_def$CHR[which(H3K27ac_count_raw_def$CHR == '1')] = 'chr1'
# H3K27ac_count_raw_def$CHR[which(H3K27ac_count_raw_def$CHR == '2')] = 'chr10'
# H3K27ac_count_raw_def$CHR[which(H3K27ac_count_raw_def$CHR == '3')] = 'chr11'
# H3K27ac_count_raw_def$CHR[which(H3K27ac_count_raw_def$CHR == '4')] = 'chr12'
# H3K27ac_count_raw_def$CHR[which(H3K27ac_count_raw_def$CHR == '5')] = 'chr13'
# H3K27ac_count_raw_def$CHR[which(H3K27ac_count_raw_def$CHR == '6')] = 'chr14'
# H3K27ac_count_raw_def$CHR[which(H3K27ac_count_raw_def$CHR == '7')] = 'chr15'
# H3K27ac_count_raw_def$CHR[which(H3K27ac_count_raw_def$CHR == '8')] = 'chr16'
# H3K27ac_count_raw_def$CHR[which(H3K27ac_count_raw_def$CHR == '9')] = 'chr17'
# H3K27ac_count_raw_def$CHR[which(H3K27ac_count_raw_def$CHR == '10')] = 'chr18'
# H3K27ac_count_raw_def$CHR[which(H3K27ac_count_raw_def$CHR == '11')] = 'chr19'
# H3K27ac_count_raw_def$CHR[which(H3K27ac_count_raw_def$CHR == '12')] = 'chr2'
# H3K27ac_count_raw_def$CHR[which(H3K27ac_count_raw_def$CHR == '13')] = 'chr3'
# H3K27ac_count_raw_def$CHR[which(H3K27ac_count_raw_def$CHR == '14')] = 'chr4'
# H3K27ac_count_raw_def$CHR[which(H3K27ac_count_raw_def$CHR == '15')] = 'chr5'
# H3K27ac_count_raw_def$CHR[which(H3K27ac_count_raw_def$CHR == '16')] = 'chr6'
# H3K27ac_count_raw_def$CHR[which(H3K27ac_count_raw_def$CHR == '17')] = 'chr7'
# H3K27ac_count_raw_def$CHR[which(H3K27ac_count_raw_def$CHR == '18')] = 'chr8'
# H3K27ac_count_raw_def$CHR[which(H3K27ac_count_raw_def$CHR == '19')] = 'chr9'
# H3K27ac_count_raw_def$CHR[which(H3K27ac_count_raw_def$CHR == '20')] = 'chrX'
# H3K27ac_count_raw_def$CHR[which(H3K27ac_count_raw_def$CHR == '21')] = 'chrY'
# chr_list = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
#              "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
#              "chrX","chrY")
# miss <- c()
# for (i in 1:nrow(H3K27ac_count_raw_def)) {
#   if(!(H3K27ac_count_raw_def$CHR[i] %in% chr_list)) miss <- append(miss,i)
# }
# # new_H3K27ac_count_def = H3K27ac_count_def[-miss,]
# H3K27ac_count_raw_def[,4:ncol(H3K27ac_count_raw_def)] = round(H3K27ac_count_raw_def[,4:ncol(H3K27ac_count_raw_def)])
# new_H3K27ac_count_raw_def = data.frame(H3K27ac_count_raw_def[,1:3],
#                                    rownames(H3K27ac_count_raw_def),
#                                    H3K27ac_count_raw_def[,4:ncol(H3K27ac_count_raw_def)])
# colnames(new_H3K27ac_count_raw_def)[1:4] = c('chr','start','end','peak_num')
# setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/count_files")
# write.csv(new_H3K27ac_count_raw_def,'H3K27ac_count_abcd_raw.csv',row.names = FALSE)
# write.table(new_H3K27ac_count_raw_def[,1:4],'H3K27ac_peak_abcd_raw.bed',sep = "\t",
#             quote = F,
#             col.names = F,
#             row.names = F)
# 
# 
# 
