library(ChIPpeakAnno)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(stringr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ChIPseeker)
rm(list = ls())
options(stringsAsFactors = F)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#######################Promoter_Peak_Annotation#######################################################
#######################BA##################################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/BA1")
Differential_BA1 = read.csv('Differential_BA1.csv')
BA1<- Differential_BA1[,c(2,3,4,5)]
write.table(BA1,'BA1.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
bedPeaksFile = "BA1.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
write.csv(peakAnno_df,'BA1_annotation.csv',row.names = FALSE)


setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/BA2")
Differential_BA2 = read.csv('Differential_BA2.csv')
BA2<- Differential_BA2[,c(2,3,4,5)]
write.table(BA2,'BA2.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
bedPeaksFile = "BA2.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
write.csv(peakAnno_df,'BA2_annotation.csv',row.names = FALSE)

setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/BA3")
Differential_BA3 = read.csv('Differential_BA3.csv')
BA3<- Differential_BA3[,c(2,3,4,5)]
write.table(BA3,'BA3.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
bedPeaksFile = "BA3.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
write.csv(peakAnno_df,'BA3_annotation.csv',row.names = FALSE)

#######################DC##################################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/D1C")
Differential_D1C = read.csv('Differential_D1C.csv')
D1C<- Differential_D1C[,c(2,3,4,5)]
write.table(D1C,'D1C.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
bedPeaksFile = "D1C.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
write.csv(peakAnno_df,'D1C_annotation.csv',row.names = FALSE)


setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/D2C")
Differential_D2C = read.csv('Differential_D2C.csv')
D2C<- Differential_D2C[,c(2,3,4,5)]
write.table(D2C,'D2C.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
bedPeaksFile = "D2C.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
write.csv(peakAnno_df,'D2C_annotation.csv',row.names = FALSE)

setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/D3C")
Differential_D3C = read.csv('Differential_D3C.csv')
D3C<- Differential_D3C[,c(2,3,4,5)]
write.table(D3C,'D3C.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
bedPeaksFile = "D3C.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
write.csv(peakAnno_df,'D3C_annotation.csv',row.names = FALSE)

#######################CA##################################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/CA1")
Differential_CA1 = read.csv('Differential_CA1.csv')
CA1<- Differential_CA1[,c(2,3,4,5)]
write.table(CA1,'CA1.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
bedPeaksFile = "CA1.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
write.csv(peakAnno_df,'CA1_annotation.csv',row.names = FALSE)

setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/CA2")
Differential_CA2 = read.csv('Differential_CA2.csv')
CA2<- Differential_CA2[,c(2,3,4,5)]
write.table(CA2,'CA2.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
bedPeaksFile = "CA2.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
write.csv(peakAnno_df,'CA2_annotation.csv',row.names = FALSE)

setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/CA3")
Differential_CA3 = read.csv('Differential_CA3.csv')
CA3<- Differential_CA3[,c(2,3,4,5)]
write.table(CA3,'CA3.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
bedPeaksFile = "CA3.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
write.csv(peakAnno_df,'CA3_annotation.csv',row.names = FALSE)

#######################BD##################################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/BD1")
Differential_BD1 = read.csv('Differential_BD1.csv')
BD1<- Differential_BD1[,c(2,3,4,5)]
write.table(BD1,'BD1.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
bedPeaksFile = "BD1.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
write.csv(peakAnno_df,'BD1_annotation.csv',row.names = FALSE)

setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/BD2")
Differential_BD2 = read.csv('Differential_BD2.csv')
BD2<- Differential_BD2[,c(2,3,4,5)]
write.table(BD2,'BD2.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
bedPeaksFile = "BD2.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
write.csv(peakAnno_df,'BD2_annotation.csv',row.names = FALSE)

setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/BD3")
Differential_BD3 = read.csv('Differential_BD3.csv')
BD3<- Differential_BD3[,c(2,3,4,5)]
write.table(BD3,'BD3.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
bedPeaksFile = "BD3.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
write.csv(peakAnno_df,'BD3_annotation.csv',row.names = FALSE)

#######################EG##################################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/EG")
Differential_EG = read.csv('Differential_EG.csv')
EG<- Differential_EG[,c(2,3,4,5)]
write.table(EG,'EG.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
bedPeaksFile = "EG.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
write.csv(peakAnno_df,'EG_annotation.csv',row.names = FALSE)

#######################EF##################################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/FE")
Differential_FE = read.csv('Differential_FE.csv')
FE<- Differential_FE[,c(2,3,4,5)]
write.table(FE,'FE.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
bedPeaksFile = "FE.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
write.csv(peakAnno_df,'FE_annotation.csv',row.names = FALSE)

#######################FG##################################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/promoter/FG")
Differential_FG = read.csv('Differential_FG.csv')
FG<- Differential_FG[,c(2,3,4,5)]
write.table(FG,'FG.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
bedPeaksFile = "FG.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
write.csv(peakAnno_df,'FG_annotation.csv',row.names = FALSE)

