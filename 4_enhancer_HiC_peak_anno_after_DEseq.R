############preparation##############
library(biomaRt)
library(dplyr)
library(DESeq2)
library(umap)
library(ggplot2)
library(ChIPpeakAnno)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(stringr)
library(reshape2)
library(ggpubr)
library(ChIPseeker)
options(stringsAsFactors = F)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#####
#######################BA1##################################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/BA1")
Differential_BA1 = read.csv('Differential_BA1.csv')
BA1<- Differential_BA1[,c(2,3,4,5)]
write.table(BA1,'BA1.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
file1 <- read.table("BA1.bed", header = FALSE)
file2 <- read.csv("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/HiC/Cortex_HiC.csv", header = FALSE)

# Create an empty data frame to store matching rows
matched_rows <- data.frame()

# Create an empty data frame to store non-matching rows from file1
non_matched_rows <- data.frame()

# Iterate through each row in file1
for (i in 1:nrow(file1)) {
  # Get the current row from file1
  row1 <- file1[i, ]
  
  # Find matching rows in file2 based on column 1
  matching_rows_file2 <- file2[file2$V1 == row1$V1, ]
  
  # Initialize a flag to determine if the row1 satisfies the conditions
  satisfied_conditions <- FALSE
  
  # Iterate through matching rows in file2
  for (j in 1:nrow(matching_rows_file2)) {
    # Get the current matching row from file2
    row2 <- matching_rows_file2[j, ]
    
    # Check conditions
    if (row1$V2 <= row2$V3 && row1$V3 >= row2$V2) {
      # Set the flag to TRUE if conditions are met
      satisfied_conditions <- TRUE
      # Add the row from file1 to the data frame if conditions are met
      # Also add column 4 of file 2 to the new data frame
      matched_row <- c(row1, row2[4])  # Combine columns from both rows
      matched_rows <- rbind(matched_rows, matched_row)
      break  # Break out of the inner loop once a match is found
    }
  }
  
  # If conditions are not met, add the row to non_matched_rows
  if (!satisfied_conditions) {
    non_matched_rows <- rbind(non_matched_rows, row1)
  }
}
# # Write the matched_rows data frame to a BED file
# write.table(matched_rows, "matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
# Write the non_matched_rows data frame to a BED file
write.table(non_matched_rows, "non_matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

bedPeaksFile = "non_matched_rows.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
unmatched_column <- peakAnno_df[, 17]
matched_column <- matched_rows[, 5]

combined_column <- c(unmatched_column, matched_column)

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
unique_genes <- unique(combined_column)
gene_map <- getBM(
  filters = "mgi_symbol",          # Adjust if your gene names are not MGI symbols
  attributes = c("mgi_symbol", "ensembl_gene_id"),
  values = unique_genes,
  mart = mart
)
gene_map_vector <- setNames(gene_map$ensembl_gene_id, gene_map$mgi_symbol)
combined_column_converted <- gene_map_vector[combined_column]
combined_data <- data.frame(Value = combined_column_converted)
write.csv(combined_data, "BA1_annotation.csv", row.names = FALSE)

#######################BA2##################################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/BA2")
Differential_BA2 = read.csv('Differential_BA2.csv')
BA2<- Differential_BA2[,c(2,3,4,5)]
write.table(BA2,'BA2.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
file1 <- read.table("BA2.bed", header = FALSE)
file2 <- read.csv("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/HiC/Cortex_HiC.csv", header = FALSE)

# Create an empty data frame to store matching rows
matched_rows <- data.frame()

# Create an empty data frame to store non-matching rows from file1
non_matched_rows <- data.frame()

# Iterate through each row in file1
for (i in 1:nrow(file1)) {
  # Get the current row from file1
  row1 <- file1[i, ]
  
  # Find matching rows in file2 based on column 1
  matching_rows_file2 <- file2[file2$V1 == row1$V1, ]
  
  # Initialize a flag to determine if the row1 satisfies the conditions
  satisfied_conditions <- FALSE
  
  # Iterate through matching rows in file2
  for (j in 1:nrow(matching_rows_file2)) {
    # Get the current matching row from file2
    row2 <- matching_rows_file2[j, ]
    
    # Check conditions
    if (row1$V2 <= row2$V3 && row1$V3 >= row2$V2) {
      # Set the flag to TRUE if conditions are met
      satisfied_conditions <- TRUE
      # Add the row from file1 to the data frame if conditions are met
      # Also add column 4 of file 2 to the new data frame
      matched_row <- c(row1, row2[4])  # Combine columns from both rows
      matched_rows <- rbind(matched_rows, matched_row)
      break  # Break out of the inner loop once a match is found
    }
  }
  
  # If conditions are not met, add the row to non_matched_rows
  if (!satisfied_conditions) {
    non_matched_rows <- rbind(non_matched_rows, row1)
  }
}
# # Write the matched_rows data frame to a BED file
# write.table(matched_rows, "matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
# Write the non_matched_rows data frame to a BED file
write.table(non_matched_rows, "non_matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

bedPeaksFile = "non_matched_rows.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
unmatched_column <- peakAnno_df[, 17]
matched_column <- matched_rows[, 5]

combined_column <- c(unmatched_column, matched_column)

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
unique_genes <- unique(combined_column)
gene_map <- getBM(
  filters = "mgi_symbol",          # Adjust if your gene names are not MGI symbols
  attributes = c("mgi_symbol", "ensembl_gene_id"),
  values = unique_genes,
  mart = mart
)
gene_map_vector <- setNames(gene_map$ensembl_gene_id, gene_map$mgi_symbol)
combined_column_converted <- gene_map_vector[combined_column]
combined_data <- data.frame(Value = combined_column_converted)
write.csv(combined_data, "BA2_annotation.csv", row.names = FALSE)

#######################BA3##################################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/BA3")
Differential_BA3 = read.csv('Differential_BA3.csv')
BA3<- Differential_BA3[,c(2,3,4,5)]
write.table(BA3,'BA3.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
file1 <- read.table("BA3.bed", header = FALSE)
file2 <- read.csv("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/HiC/Cortex_HiC.csv", header = FALSE)

# Create an empty data frame to store matching rows
matched_rows <- data.frame()

# Create an empty data frame to store non-matching rows from file1
non_matched_rows <- data.frame()

# Iterate through each row in file1
for (i in 1:nrow(file1)) {
  # Get the current row from file1
  row1 <- file1[i, ]
  
  # Find matching rows in file2 based on column 1
  matching_rows_file2 <- file2[file2$V1 == row1$V1, ]
  
  # Initialize a flag to determine if the row1 satisfies the conditions
  satisfied_conditions <- FALSE
  
  # Iterate through matching rows in file2
  for (j in 1:nrow(matching_rows_file2)) {
    # Get the current matching row from file2
    row2 <- matching_rows_file2[j, ]
    
    # Check conditions
    if (row1$V2 <= row2$V3 && row1$V3 >= row2$V2) {
      # Set the flag to TRUE if conditions are met
      satisfied_conditions <- TRUE
      # Add the row from file1 to the data frame if conditions are met
      # Also add column 4 of file 2 to the new data frame
      matched_row <- c(row1, row2[4])  # Combine columns from both rows
      matched_rows <- rbind(matched_rows, matched_row)
      break  # Break out of the inner loop once a match is found
    }
  }
  
  # If conditions are not met, add the row to non_matched_rows
  if (!satisfied_conditions) {
    non_matched_rows <- rbind(non_matched_rows, row1)
  }
}
# # Write the matched_rows data frame to a BED file
# write.table(matched_rows, "matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
# Write the non_matched_rows data frame to a BED file
write.table(non_matched_rows, "non_matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

bedPeaksFile = "non_matched_rows.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
unmatched_column <- peakAnno_df[, 17]
matched_column <- matched_rows[, 5]

combined_column <- c(unmatched_column, matched_column)

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
unique_genes <- unique(combined_column)
gene_map <- getBM(
  filters = "mgi_symbol",          # Adjust if your gene names are not MGI symbols
  attributes = c("mgi_symbol", "ensembl_gene_id"),
  values = unique_genes,
  mart = mart
)
gene_map_vector <- setNames(gene_map$ensembl_gene_id, gene_map$mgi_symbol)
combined_column_converted <- gene_map_vector[combined_column]
combined_data <- data.frame(Value = combined_column_converted)
write.csv(combined_data, "BA3_annotation.csv", row.names = FALSE)
#######################D1C##################################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/D1C")
Differential_D1C = read.csv('Differential_D1C.csv')
D1C<- Differential_D1C[,c(2,3,4,5)]
write.table(D1C,'D1C.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
file1 <- read.table("D1C.bed", header = FALSE)
file2 <- read.csv("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/HiC/Cortex_HiC.csv", header = FALSE)

# Create an empty data frame to store matching rows
matched_rows <- data.frame()

# Create an empty data frame to store non-matching rows from file1
non_matched_rows <- data.frame()

# Iterate through each row in file1
for (i in 1:nrow(file1)) {
  # Get the current row from file1
  row1 <- file1[i, ]
  
  # Find matching rows in file2 based on column 1
  matching_rows_file2 <- file2[file2$V1 == row1$V1, ]
  
  # Initialize a flag to determine if the row1 satisfies the conditions
  satisfied_conditions <- FALSE
  
  # Iterate through matching rows in file2
  for (j in 1:nrow(matching_rows_file2)) {
    # Get the current matching row from file2
    row2 <- matching_rows_file2[j, ]
    
    # Check conditions
    if (row1$V2 <= row2$V3 && row1$V3 >= row2$V2) {
      # Set the flag to TRUE if conditions are met
      satisfied_conditions <- TRUE
      # Add the row from file1 to the data frame if conditions are met
      # Also add column 4 of file 2 to the new data frame
      matched_row <- c(row1, row2[4])  # Combine columns from both rows
      matched_rows <- rbind(matched_rows, matched_row)
      break  # Break out of the inner loop once a match is found
    }
  }
  
  # If conditions are not met, add the row to non_matched_rows
  if (!satisfied_conditions) {
    non_matched_rows <- rbind(non_matched_rows, row1)
  }
}
# # Write the matched_rows data frame to a BED file
# write.table(matched_rows, "matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
# Write the non_matched_rows data frame to a BED file
write.table(non_matched_rows, "non_matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

bedPeaksFile = "non_matched_rows.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
unmatched_column <- peakAnno_df[, 17]
matched_column <- matched_rows[, 5]

combined_column <- c(unmatched_column, matched_column)

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
unique_genes <- unique(combined_column)
gene_map <- getBM(
  filters = "mgi_symbol",          # Adjust if your gene names are not MGI symbols
  attributes = c("mgi_symbol", "ensembl_gene_id"),
  values = unique_genes,
  mart = mart
)
gene_map_vector <- setNames(gene_map$ensembl_gene_id, gene_map$mgi_symbol)
combined_column_converted <- gene_map_vector[combined_column]
combined_data <- data.frame(Value = combined_column_converted)
write.csv(combined_data, "D1C_annotation.csv", row.names = FALSE)

#######################D2C##################################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/D2C")
Differential_D2C = read.csv('Differential_D2C.csv')
D2C<- Differential_D2C[,c(2,3,4,5)]
write.table(D2C,'D2C.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
file1 <- read.table("D2C.bed", header = FALSE)
file2 <- read.csv("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/HiC/Cortex_HiC.csv", header = FALSE)

# Create an empty data frame to store matching rows
matched_rows <- data.frame()

# Create an empty data frame to store non-matching rows from file1
non_matched_rows <- data.frame()

# Iterate through each row in file1
for (i in 1:nrow(file1)) {
  # Get the current row from file1
  row1 <- file1[i, ]
  
  # Find matching rows in file2 based on column 1
  matching_rows_file2 <- file2[file2$V1 == row1$V1, ]
  
  # Initialize a flag to determine if the row1 satisfies the conditions
  satisfied_conditions <- FALSE
  
  # Iterate through matching rows in file2
  for (j in 1:nrow(matching_rows_file2)) {
    # Get the current matching row from file2
    row2 <- matching_rows_file2[j, ]
    
    # Check conditions
    if (row1$V2 <= row2$V3 && row1$V3 >= row2$V2) {
      # Set the flag to TRUE if conditions are met
      satisfied_conditions <- TRUE
      # Add the row from file1 to the data frame if conditions are met
      # Also add column 4 of file 2 to the new data frame
      matched_row <- c(row1, row2[4])  # Combine columns from both rows
      matched_rows <- rbind(matched_rows, matched_row)
      break  # Break out of the inner loop once a match is found
    }
  }
  
  # If conditions are not met, add the row to non_matched_rows
  if (!satisfied_conditions) {
    non_matched_rows <- rbind(non_matched_rows, row1)
  }
}
# # Write the matched_rows data frame to a BED file
# write.table(matched_rows, "matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
# Write the non_matched_rows data frame to a BED file
write.table(non_matched_rows, "non_matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

bedPeaksFile = "non_matched_rows.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
unmatched_column <- peakAnno_df[, 17]
matched_column <- matched_rows[, 5]

combined_column <- c(unmatched_column, matched_column)

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
unique_genes <- unique(combined_column)
gene_map <- getBM(
  filters = "mgi_symbol",          # Adjust if your gene names are not MGI symbols
  attributes = c("mgi_symbol", "ensembl_gene_id"),
  values = unique_genes,
  mart = mart
)
gene_map_vector <- setNames(gene_map$ensembl_gene_id, gene_map$mgi_symbol)
combined_column_converted <- gene_map_vector[combined_column]
combined_data <- data.frame(Value = combined_column_converted)
write.csv(combined_data, "D2C_annotation.csv", row.names = FALSE)

#######################D3C##################################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/D3C")
Differential_D3C = read.csv('Differential_D3C.csv')
D3C<- Differential_D3C[,c(2,3,4,5)]
write.table(D3C,'D3C.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
file1 <- read.table("D3C.bed", header = FALSE)
file2 <- read.csv("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/HiC/Cortex_HiC.csv", header = FALSE)

# Create an empty data frame to store matching rows
matched_rows <- data.frame()

# Create an empty data frame to store non-matching rows from file1
non_matched_rows <- data.frame()

# Iterate through each row in file1
for (i in 1:nrow(file1)) {
  # Get the current row from file1
  row1 <- file1[i, ]
  
  # Find matching rows in file2 based on column 1
  matching_rows_file2 <- file2[file2$V1 == row1$V1, ]
  
  # Initialize a flag to determine if the row1 satisfies the conditions
  satisfied_conditions <- FALSE
  
  # Iterate through matching rows in file2
  for (j in 1:nrow(matching_rows_file2)) {
    # Get the current matching row from file2
    row2 <- matching_rows_file2[j, ]
    
    # Check conditions
    if (row1$V2 <= row2$V3 && row1$V3 >= row2$V2) {
      # Set the flag to TRUE if conditions are met
      satisfied_conditions <- TRUE
      # Add the row from file1 to the data frame if conditions are met
      # Also add column 4 of file 2 to the new data frame
      matched_row <- c(row1, row2[4])  # Combine columns from both rows
      matched_rows <- rbind(matched_rows, matched_row)
      break  # Break out of the inner loop once a match is found
    }
  }
  
  # If conditions are not met, add the row to non_matched_rows
  if (!satisfied_conditions) {
    non_matched_rows <- rbind(non_matched_rows, row1)
  }
}
# # Write the matched_rows data frame to a BED file
# write.table(matched_rows, "matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
# Write the non_matched_rows data frame to a BED file
write.table(non_matched_rows, "non_matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

bedPeaksFile = "non_matched_rows.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
unmatched_column <- peakAnno_df[, 17]
matched_column <- matched_rows[, 5]

combined_column <- c(unmatched_column, matched_column)

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
unique_genes <- unique(combined_column)
gene_map <- getBM(
  filters = "mgi_symbol",          # Adjust if your gene names are not MGI symbols
  attributes = c("mgi_symbol", "ensembl_gene_id"),
  values = unique_genes,
  mart = mart
)
gene_map_vector <- setNames(gene_map$ensembl_gene_id, gene_map$mgi_symbol)
combined_column_converted <- gene_map_vector[combined_column]
combined_data <- data.frame(Value = combined_column_converted)
write.csv(combined_data, "D3C_annotation.csv", row.names = FALSE)
#######################CA1##################################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/CA1")
Differential_CA1 = read.csv('Differential_CA1.csv')
CA1<- Differential_CA1[,c(2,3,4,5)]
write.table(CA1,'CA1.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
file1 <- read.table("CA1.bed", header = FALSE)
file2 <- read.csv("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/HiC/Cortex_HiC.csv", header = FALSE)

# Create an empty data frame to store matching rows
matched_rows <- data.frame()

# Create an empty data frame to store non-matching rows from file1
non_matched_rows <- data.frame()

# Iterate through each row in file1
for (i in 1:nrow(file1)) {
  # Get the current row from file1
  row1 <- file1[i, ]
  
  # Find matching rows in file2 CAsed on column 1
  matching_rows_file2 <- file2[file2$V1 == row1$V1, ]
  
  # Initialize a flag to determine if the row1 satisfies the conditions
  satisfied_conditions <- FALSE
  
  # Iterate through matching rows in file2
  for (j in 1:nrow(matching_rows_file2)) {
    # Get the current matching row from file2
    row2 <- matching_rows_file2[j, ]
    
    # Check conditions
    if (row1$V2 <= row2$V3 && row1$V3 >= row2$V2) {
      # Set the flag to TRUE if conditions are met
      satisfied_conditions <- TRUE
      # Add the row from file1 to the data frame if conditions are met
      # Also add column 4 of file 2 to the new data frame
      matched_row <- c(row1, row2[4])  # Combine columns from both rows
      matched_rows <- rbind(matched_rows, matched_row)
      break  # Break out of the inner loop once a match is found
    }
  }
  
  # If conditions are not met, add the row to non_matched_rows
  if (!satisfied_conditions) {
    non_matched_rows <- rbind(non_matched_rows, row1)
  }
}
# # Write the matched_rows data frame to a BED file
# write.table(matched_rows, "matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
# Write the non_matched_rows data frame to a BED file
write.table(non_matched_rows, "non_matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

bedPeaksFile = "non_matched_rows.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
unmatched_column <- peakAnno_df[, 17]
matched_column <- matched_rows[, 5]

combined_column <- c(unmatched_column, matched_column)

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
unique_genes <- unique(combined_column)
gene_map <- getBM(
  filters = "mgi_symbol",          # Adjust if your gene names are not MGI symbols
  attributes = c("mgi_symbol", "ensembl_gene_id"),
  values = unique_genes,
  mart = mart
)
gene_map_vector <- setNames(gene_map$ensembl_gene_id, gene_map$mgi_symbol)
combined_column_converted <- gene_map_vector[combined_column]
combined_data <- data.frame(Value = combined_column_converted)
write.csv(combined_data, "CA1_annotation.csv", row.names = FALSE)

#######################CA2##################################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/CA2")
Differential_CA2 = read.csv('Differential_CA2.csv')
CA2<- Differential_CA2[,c(2,3,4,5)]
write.table(CA2,'CA2.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
file1 <- read.table("CA2.bed", header = FALSE)
file2 <- read.csv("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/HiC/Cortex_HiC.csv", header = FALSE)

# Create an empty data frame to store matching rows
matched_rows <- data.frame()

# Create an empty data frame to store non-matching rows from file1
non_matched_rows <- data.frame()

# Iterate through each row in file1
for (i in 1:nrow(file1)) {
  # Get the current row from file1
  row1 <- file1[i, ]
  
  # Find matching rows in file2 CAsed on column 1
  matching_rows_file2 <- file2[file2$V1 == row1$V1, ]
  
  # Initialize a flag to determine if the row1 satisfies the conditions
  satisfied_conditions <- FALSE
  
  # Iterate through matching rows in file2
  for (j in 1:nrow(matching_rows_file2)) {
    # Get the current matching row from file2
    row2 <- matching_rows_file2[j, ]
    
    # Check conditions
    if (row1$V2 <= row2$V3 && row1$V3 >= row2$V2) {
      # Set the flag to TRUE if conditions are met
      satisfied_conditions <- TRUE
      # Add the row from file1 to the data frame if conditions are met
      # Also add column 4 of file 2 to the new data frame
      matched_row <- c(row1, row2[4])  # Combine columns from both rows
      matched_rows <- rbind(matched_rows, matched_row)
      break  # Break out of the inner loop once a match is found
    }
  }
  
  # If conditions are not met, add the row to non_matched_rows
  if (!satisfied_conditions) {
    non_matched_rows <- rbind(non_matched_rows, row1)
  }
}
# # Write the matched_rows data frame to a BED file
# write.table(matched_rows, "matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
# Write the non_matched_rows data frame to a BED file
write.table(non_matched_rows, "non_matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

bedPeaksFile = "non_matched_rows.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
unmatched_column <- peakAnno_df[, 17]
matched_column <- matched_rows[, 5]

combined_column <- c(unmatched_column, matched_column)

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
unique_genes <- unique(combined_column)
gene_map <- getBM(
  filters = "mgi_symbol",          # Adjust if your gene names are not MGI symbols
  attributes = c("mgi_symbol", "ensembl_gene_id"),
  values = unique_genes,
  mart = mart
)
gene_map_vector <- setNames(gene_map$ensembl_gene_id, gene_map$mgi_symbol)
combined_column_converted <- gene_map_vector[combined_column]
combined_data <- data.frame(Value = combined_column_converted)
write.csv(combined_data, "CA2_annotation.csv", row.names = FALSE)

#######################CA3##################################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/CA3")
Differential_CA3 = read.csv('Differential_CA3.csv')
CA3<- Differential_CA3[,c(2,3,4,5)]
write.table(CA3,'CA3.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
file1 <- read.table("CA3.bed", header = FALSE)
file2 <- read.csv("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/HiC/Cortex_HiC.csv", header = FALSE)

# Create an empty data frame to store matching rows
matched_rows <- data.frame()

# Create an empty data frame to store non-matching rows from file1
non_matched_rows <- data.frame()

# Iterate through each row in file1
for (i in 1:nrow(file1)) {
  # Get the current row from file1
  row1 <- file1[i, ]
  
  # Find matching rows in file2 CAsed on column 1
  matching_rows_file2 <- file2[file2$V1 == row1$V1, ]
  
  # Initialize a flag to determine if the row1 satisfies the conditions
  satisfied_conditions <- FALSE
  
  # Iterate through matching rows in file2
  for (j in 1:nrow(matching_rows_file2)) {
    # Get the current matching row from file2
    row2 <- matching_rows_file2[j, ]
    
    # Check conditions
    if (row1$V2 <= row2$V3 && row1$V3 >= row2$V2) {
      # Set the flag to TRUE if conditions are met
      satisfied_conditions <- TRUE
      # Add the row from file1 to the data frame if conditions are met
      # Also add column 4 of file 2 to the new data frame
      matched_row <- c(row1, row2[4])  # Combine columns from both rows
      matched_rows <- rbind(matched_rows, matched_row)
      break  # Break out of the inner loop once a match is found
    }
  }
  
  # If conditions are not met, add the row to non_matched_rows
  if (!satisfied_conditions) {
    non_matched_rows <- rbind(non_matched_rows, row1)
  }
}
# # Write the matched_rows data frame to a BED file
# write.table(matched_rows, "matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
# Write the non_matched_rows data frame to a BED file
write.table(non_matched_rows, "non_matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

bedPeaksFile = "non_matched_rows.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
unmatched_column <- peakAnno_df[, 17]
matched_column <- matched_rows[, 5]

combined_column <- c(unmatched_column, matched_column)

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
unique_genes <- unique(combined_column)
gene_map <- getBM(
  filters = "mgi_symbol",          # Adjust if your gene names are not MGI symbols
  attributes = c("mgi_symbol", "ensembl_gene_id"),
  values = unique_genes,
  mart = mart
)
gene_map_vector <- setNames(gene_map$ensembl_gene_id, gene_map$mgi_symbol)
combined_column_converted <- gene_map_vector[combined_column]
combined_data <- data.frame(Value = combined_column_converted)
write.csv(combined_data, "CA3_annotation.csv", row.names = FALSE)
#######################BD1##################################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/BD1")
Differential_BD1 = read.csv('Differential_BD1.csv')
BD1<- Differential_BD1[,c(2,3,4,5)]
write.table(BD1,'BD1.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
file1 <- read.table("BD1.bed", header = FALSE)
file2 <- read.csv("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/HiC/Cortex_HiC.csv", header = FALSE)

# Create an empty data frame to store matching rows
matched_rows <- data.frame()

# Create an empty data frame to store non-matching rows from file1
non_matched_rows <- data.frame()

# Iterate through each row in file1
for (i in 1:nrow(file1)) {
  # Get the current row from file1
  row1 <- file1[i, ]
  
  # Find matching rows in file2 BDsed on column 1
  matching_rows_file2 <- file2[file2$V1 == row1$V1, ]
  
  # Initialize a flag to determine if the row1 satisfies the conditions
  satisfied_conditions <- FALSE
  
  # Iterate through matching rows in file2
  for (j in 1:nrow(matching_rows_file2)) {
    # Get the current matching row from file2
    row2 <- matching_rows_file2[j, ]
    
    # Check conditions
    if (row1$V2 <= row2$V3 && row1$V3 >= row2$V2) {
      # Set the flag to TRUE if conditions are met
      satisfied_conditions <- TRUE
      # Add the row from file1 to the data frame if conditions are met
      # Also add column 4 of file 2 to the new data frame
      matched_row <- c(row1, row2[4])  # Combine columns from both rows
      matched_rows <- rbind(matched_rows, matched_row)
      break  # Break out of the inner loop once a match is found
    }
  }
  
  # If conditions are not met, add the row to non_matched_rows
  if (!satisfied_conditions) {
    non_matched_rows <- rbind(non_matched_rows, row1)
  }
}
# # Write the matched_rows data frame to a BED file
# write.table(matched_rows, "matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
# Write the non_matched_rows data frame to a BED file
write.table(non_matched_rows, "non_matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

bedPeaksFile = "non_matched_rows.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
unmatched_column <- peakAnno_df[, 17]
matched_column <- matched_rows[, 5]

combined_column <- c(unmatched_column, matched_column)

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
unique_genes <- unique(combined_column)
gene_map <- getBM(
  filters = "mgi_symbol",          # Adjust if your gene names are not MGI symbols
  attributes = c("mgi_symbol", "ensembl_gene_id"),
  values = unique_genes,
  mart = mart
)
gene_map_vector <- setNames(gene_map$ensembl_gene_id, gene_map$mgi_symbol)
combined_column_converted <- gene_map_vector[combined_column]
combined_data <- data.frame(Value = combined_column_converted)
write.csv(combined_data, "BD1_annotation.csv", row.names = FALSE)

#######################BD2##################################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/BD2")
Differential_BD2 = read.csv('Differential_BD2.csv')
BD2<- Differential_BD2[,c(2,3,4,5)]
write.table(BD2,'BD2.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
file1 <- read.table("BD2.bed", header = FALSE)
file2 <- read.csv("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/HiC/Cortex_HiC.csv", header = FALSE)

# Create an empty data frame to store matching rows
matched_rows <- data.frame()

# Create an empty data frame to store non-matching rows from file1
non_matched_rows <- data.frame()

# Iterate through each row in file1
for (i in 1:nrow(file1)) {
  # Get the current row from file1
  row1 <- file1[i, ]
  
  # Find matching rows in file2 BDsed on column 1
  matching_rows_file2 <- file2[file2$V1 == row1$V1, ]
  
  # Initialize a flag to determine if the row1 satisfies the conditions
  satisfied_conditions <- FALSE
  
  # Iterate through matching rows in file2
  for (j in 1:nrow(matching_rows_file2)) {
    # Get the current matching row from file2
    row2 <- matching_rows_file2[j, ]
    
    # Check conditions
    if (row1$V2 <= row2$V3 && row1$V3 >= row2$V2) {
      # Set the flag to TRUE if conditions are met
      satisfied_conditions <- TRUE
      # Add the row from file1 to the data frame if conditions are met
      # Also add column 4 of file 2 to the new data frame
      matched_row <- c(row1, row2[4])  # Combine columns from both rows
      matched_rows <- rbind(matched_rows, matched_row)
      break  # Break out of the inner loop once a match is found
    }
  }
  
  # If conditions are not met, add the row to non_matched_rows
  if (!satisfied_conditions) {
    non_matched_rows <- rbind(non_matched_rows, row1)
  }
}
# # Write the matched_rows data frame to a BED file
# write.table(matched_rows, "matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
# Write the non_matched_rows data frame to a BED file
write.table(non_matched_rows, "non_matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

bedPeaksFile = "non_matched_rows.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
unmatched_column <- peakAnno_df[, 17]
matched_column <- matched_rows[, 5]

combined_column <- c(unmatched_column, matched_column)

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
unique_genes <- unique(combined_column)
gene_map <- getBM(
  filters = "mgi_symbol",          # Adjust if your gene names are not MGI symbols
  attributes = c("mgi_symbol", "ensembl_gene_id"),
  values = unique_genes,
  mart = mart
)
gene_map_vector <- setNames(gene_map$ensembl_gene_id, gene_map$mgi_symbol)
combined_column_converted <- gene_map_vector[combined_column]
combined_data <- data.frame(Value = combined_column_converted)
write.csv(combined_data, "BD2_annotation.csv", row.names = FALSE)

#######################BD3##################################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/BD3")
Differential_BD3 = read.csv('Differential_BD3.csv')
BD3<- Differential_BD3[,c(2,3,4,5)]
write.table(BD3,'BD3.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
file1 <- read.table("BD3.bed", header = FALSE)
file2 <- read.csv("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/HiC/Cortex_HiC.csv", header = FALSE)

# Create an empty data frame to store matching rows
matched_rows <- data.frame()

# Create an empty data frame to store non-matching rows from file1
non_matched_rows <- data.frame()

# Iterate through each row in file1
for (i in 1:nrow(file1)) {
  # Get the current row from file1
  row1 <- file1[i, ]
  
  # Find matching rows in file2 BDsed on column 1
  matching_rows_file2 <- file2[file2$V1 == row1$V1, ]
  
  # Initialize a flag to determine if the row1 satisfies the conditions
  satisfied_conditions <- FALSE
  
  # Iterate through matching rows in file2
  for (j in 1:nrow(matching_rows_file2)) {
    # Get the current matching row from file2
    row2 <- matching_rows_file2[j, ]
    
    # Check conditions
    if (row1$V2 <= row2$V3 && row1$V3 >= row2$V2) {
      # Set the flag to TRUE if conditions are met
      satisfied_conditions <- TRUE
      # Add the row from file1 to the data frame if conditions are met
      # Also add column 4 of file 2 to the new data frame
      matched_row <- c(row1, row2[4])  # Combine columns from both rows
      matched_rows <- rbind(matched_rows, matched_row)
      break  # Break out of the inner loop once a match is found
    }
  }
  
  # If conditions are not met, add the row to non_matched_rows
  if (!satisfied_conditions) {
    non_matched_rows <- rbind(non_matched_rows, row1)
  }
}
# # Write the matched_rows data frame to a BED file
# write.table(matched_rows, "matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
# Write the non_matched_rows data frame to a BED file
write.table(non_matched_rows, "non_matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

bedPeaksFile = "non_matched_rows.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
unmatched_column <- peakAnno_df[, 17]
matched_column <- matched_rows[, 5]

combined_column <- c(unmatched_column, matched_column)

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
unique_genes <- unique(combined_column)
gene_map <- getBM(
  filters = "mgi_symbol",          # Adjust if your gene names are not MGI symbols
  attributes = c("mgi_symbol", "ensembl_gene_id"),
  values = unique_genes,
  mart = mart
)
gene_map_vector <- setNames(gene_map$ensembl_gene_id, gene_map$mgi_symbol)
combined_column_converted <- gene_map_vector[combined_column]
combined_data <- data.frame(Value = combined_column_converted)
write.csv(combined_data, "BD3_annotation.csv", row.names = FALSE)
#######################A1E##################################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/A1E")
Differential_A1E = read.csv('Differential_A1E.csv')
A1E<- Differential_A1E[,c(2,3,4,5)]
write.table(A1E,'A1E.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
file1 <- read.table("A1E.bed", header = FALSE)
file2 <- read.csv("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/HiC/Cortex_HiC.csv", header = FALSE)

# Create an empty data frame to store matching rows
matched_rows <- data.frame()

# Create an empty data frame to store non-matching rows from file1
non_matched_rows <- data.frame()

# Iterate through each row in file1
for (i in 1:nrow(file1)) {
  # Get the current row from file1
  row1 <- file1[i, ]
  
  # Find matching rows in file2 based on column 1
  matching_rows_file2 <- file2[file2$V1 == row1$V1, ]
  
  # Initialize a flag to determine if the row1 satisfies the conditions
  satisfied_conditions <- FALSE
  
  # Iterate through matching rows in file2
  for (j in 1:nrow(matching_rows_file2)) {
    # Get the current matching row from file2
    row2 <- matching_rows_file2[j, ]
    
    # Check conditions
    if (row1$V2 <= row2$V3 && row1$V3 >= row2$V2) {
      # Set the flag to TRUE if conditions are met
      satisfied_conditions <- TRUE
      # Add the row from file1 to the data frame if conditions are met
      # Also add column 4 of file 2 to the new data frame
      matched_row <- c(row1, row2[4])  # Combine columns from both rows
      matched_rows <- rbind(matched_rows, matched_row)
      break  # Break out of the inner loop once a match is found
    }
  }
  
  # If conditions are not met, add the row to non_matched_rows
  if (!satisfied_conditions) {
    non_matched_rows <- rbind(non_matched_rows, row1)
  }
}
# # Write the matched_rows data frame to a BED file
# write.table(matched_rows, "matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
# Write the non_matched_rows data frame to a BED file
write.table(non_matched_rows, "non_matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

bedPeaksFile = "non_matched_rows.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
unmatched_column <- peakAnno_df[, 17]
matched_column <- matched_rows[, 5]

combined_column <- c(unmatched_column, matched_column)

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
unique_genes <- unique(combined_column)
gene_map <- getBM(
  filters = "mgi_symbol",          # Adjust if your gene names are not MGI symbols
  attributes = c("mgi_symbol", "ensembl_gene_id"),
  values = unique_genes,
  mart = mart
)
gene_map_vector <- setNames(gene_map$ensembl_gene_id, gene_map$mgi_symbol)
combined_column_converted <- gene_map_vector[combined_column]
combined_data <- data.frame(Value = combined_column_converted)
write.csv(combined_data, "A1E_annotation.csv", row.names = FALSE)

#######################A2E##################################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/A2E")
Differential_A2E = read.csv('Differential_A2E.csv')
A2E<- Differential_A2E[,c(2,3,4,5)]
write.table(A2E,'A2E.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
file1 <- read.table("A2E.bed", header = FALSE)
file2 <- read.csv("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/HiC/Cortex_HiC.csv", header = FALSE)

# Create an empty data frame to store matching rows
matched_rows <- data.frame()

# Create an empty data frame to store non-matching rows from file1
non_matched_rows <- data.frame()

# Iterate through each row in file1
for (i in 1:nrow(file1)) {
  # Get the current row from file1
  row1 <- file1[i, ]
  
  # Find matching rows in file2 based on column 1
  matching_rows_file2 <- file2[file2$V1 == row1$V1, ]
  
  # Initialize a flag to determine if the row1 satisfies the conditions
  satisfied_conditions <- FALSE
  
  # Iterate through matching rows in file2
  for (j in 1:nrow(matching_rows_file2)) {
    # Get the current matching row from file2
    row2 <- matching_rows_file2[j, ]
    
    # Check conditions
    if (row1$V2 <= row2$V3 && row1$V3 >= row2$V2) {
      # Set the flag to TRUE if conditions are met
      satisfied_conditions <- TRUE
      # Add the row from file1 to the data frame if conditions are met
      # Also add column 4 of file 2 to the new data frame
      matched_row <- c(row1, row2[4])  # Combine columns from both rows
      matched_rows <- rbind(matched_rows, matched_row)
      break  # Break out of the inner loop once a match is found
    }
  }
  
  # If conditions are not met, add the row to non_matched_rows
  if (!satisfied_conditions) {
    non_matched_rows <- rbind(non_matched_rows, row1)
  }
}
# # Write the matched_rows data frame to a BED file
# write.table(matched_rows, "matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
# Write the non_matched_rows data frame to a BED file
write.table(non_matched_rows, "non_matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

bedPeaksFile = "non_matched_rows.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
unmatched_column <- peakAnno_df[, 17]
matched_column <- matched_rows[, 5]

combined_column <- c(unmatched_column, matched_column)

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
unique_genes <- unique(combined_column)
gene_map <- getBM(
  filters = "mgi_symbol",          # Adjust if your gene names are not MGI symbols
  attributes = c("mgi_symbol", "ensembl_gene_id"),
  values = unique_genes,
  mart = mart
)
gene_map_vector <- setNames(gene_map$ensembl_gene_id, gene_map$mgi_symbol)
combined_column_converted <- gene_map_vector[combined_column]
combined_data <- data.frame(Value = combined_column_converted)
write.csv(combined_data, "A2E_annotation.csv", row.names = FALSE)

#######################A3E##################################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/A3E")
Differential_A3E = read.csv('Differential_A3E.csv')
A3E<- Differential_A3E[,c(2,3,4,5)]
write.table(A3E,'A3E.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
file1 <- read.table("A3E.bed", header = FALSE)
file2 <- read.csv("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/HiC/Cortex_HiC.csv", header = FALSE)

# Create an empty data frame to store matching rows
matched_rows <- data.frame()

# Create an empty data frame to store non-matching rows from file1
non_matched_rows <- data.frame()

# Iterate through each row in file1
for (i in 1:nrow(file1)) {
  # Get the current row from file1
  row1 <- file1[i, ]
  
  # Find matching rows in file2 based on column 1
  matching_rows_file2 <- file2[file2$V1 == row1$V1, ]
  
  # Initialize a flag to determine if the row1 satisfies the conditions
  satisfied_conditions <- FALSE
  
  # Iterate through matching rows in file2
  for (j in 1:nrow(matching_rows_file2)) {
    # Get the current matching row from file2
    row2 <- matching_rows_file2[j, ]
    
    # Check conditions
    if (row1$V2 <= row2$V3 && row1$V3 >= row2$V2) {
      # Set the flag to TRUE if conditions are met
      satisfied_conditions <- TRUE
      # Add the row from file1 to the data frame if conditions are met
      # Also add column 4 of file 2 to the new data frame
      matched_row <- c(row1, row2[4])  # Combine columns from both rows
      matched_rows <- rbind(matched_rows, matched_row)
      break  # Break out of the inner loop once a match is found
    }
  }
  
  # If conditions are not met, add the row to non_matched_rows
  if (!satisfied_conditions) {
    non_matched_rows <- rbind(non_matched_rows, row1)
  }
}
# # Write the matched_rows data frame to a BED file
# write.table(matched_rows, "matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
# Write the non_matched_rows data frame to a BED file
write.table(non_matched_rows, "non_matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

bedPeaksFile = "non_matched_rows.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
unmatched_column <- peakAnno_df[, 17]
matched_column <- matched_rows[, 5]

combined_column <- c(unmatched_column, matched_column)

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
unique_genes <- unique(combined_column)
gene_map <- getBM(
  filters = "mgi_symbol",          # Adjust if your gene names are not MGI symbols
  attributes = c("mgi_symbol", "ensembl_gene_id"),
  values = unique_genes,
  mart = mart
)
gene_map_vector <- setNames(gene_map$ensembl_gene_id, gene_map$mgi_symbol)
combined_column_converted <- gene_map_vector[combined_column]
combined_data <- data.frame(Value = combined_column_converted)
write.csv(combined_data, "A3E_annotation.csv", row.names = FALSE)
#######################EG##################################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/EG")
Differential_EG = read.csv('Differential_EG.csv')
EG<- Differential_EG[,c(2,3,4,5)]
write.table(EG,'EG.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
file1 <- read.table("EG.bed", header = FALSE)
file2 <- read.csv("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/HiC/Cortex_HiC.csv", header = FALSE)

# Create an empty data frame to store matching rows
matched_rows <- data.frame()

# Create an empty data frame to store non-matching rows from file1
non_matched_rows <- data.frame()

# Iterate through each row in file1
for (i in 1:nrow(file1)) {
  # Get the current row from file1
  row1 <- file1[i, ]
  
  # Find matching rows in file2 BDsed on column 1
  matching_rows_file2 <- file2[file2$V1 == row1$V1, ]
  
  # Initialize a flag to determine if the row1 satisfies the conditions
  satisfied_conditions <- FALSE
  
  # Iterate through matching rows in file2
  for (j in 1:nrow(matching_rows_file2)) {
    # Get the current matching row from file2
    row2 <- matching_rows_file2[j, ]
    
    # Check conditions
    if (row1$V2 <= row2$V3 && row1$V3 >= row2$V2) {
      # Set the flag to TRUE if conditions are met
      satisfied_conditions <- TRUE
      # Add the row from file1 to the data frame if conditions are met
      # Also add column 4 of file 2 to the new data frame
      matched_row <- c(row1, row2[4])  # Combine columns from both rows
      matched_rows <- rbind(matched_rows, matched_row)
      break  # Break out of the inner loop once a match is found
    }
  }
  
  # If conditions are not met, add the row to non_matched_rows
  if (!satisfied_conditions) {
    non_matched_rows <- rbind(non_matched_rows, row1)
  }
}
# # Write the matched_rows data frame to a BED file
# write.table(matched_rows, "matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
# Write the non_matched_rows data frame to a BED file
write.table(non_matched_rows, "non_matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

bedPeaksFile = "non_matched_rows.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
unmatched_column <- peakAnno_df[, 17]
matched_column <- matched_rows[, 5]

combined_column <- c(unmatched_column, matched_column)

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
unique_genes <- unique(combined_column)
gene_map <- getBM(
  filters = "mgi_symbol",          # Adjust if your gene names are not MGI symbols
  attributes = c("mgi_symbol", "ensembl_gene_id"),
  values = unique_genes,
  mart = mart
)
gene_map_vector <- setNames(gene_map$ensembl_gene_id, gene_map$mgi_symbol)
combined_column_converted <- gene_map_vector[combined_column]
combined_data <- data.frame(Value = combined_column_converted)
write.csv(combined_data, "EG_annotation.csv", row.names = FALSE)
#######################FE##################################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/FE")
Differential_FE = read.csv('Differential_FE.csv')
FE<- Differential_FE[,c(2,3,4,5)]
write.table(FE,'FE.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
file1 <- read.table("FE.bed", header = FALSE)
file2 <- read.csv("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/HiC/Cortex_HiC.csv", header = FALSE)

# Create an empty data frame to store matching rows
matched_rows <- data.frame()

# Create an empty data frame to store non-matching rows from file1
non_matched_rows <- data.frame()

# Iterate through each row in file1
for (i in 1:nrow(file1)) {
  # Get the current row from file1
  row1 <- file1[i, ]
  
  # Find matching rows in file2 based on column 1
  matching_rows_file2 <- file2[file2$V1 == row1$V1, ]
  
  # Initialize a flag to determine if the row1 satisfies the conditions
  satisfied_conditions <- FALSE
  
  # Iterate through matching rows in file2
  for (j in 1:nrow(matching_rows_file2)) {
    # Get the current matching row from file2
    row2 <- matching_rows_file2[j, ]
    
    # Check conditions
    if (row1$V2 <= row2$V3 && row1$V3 >= row2$V2) {
      # Set the flag to TRUE if conditions are met
      satisfied_conditions <- TRUE
      # Add the row from file1 to the data frame if conditions are met
      # Also add column 4 of file 2 to the new data frame
      matched_row <- c(row1, row2[4])  # Combine columns from both rows
      matched_rows <- rbind(matched_rows, matched_row)
      break  # Break out of the inner loop once a match is found
    }
  }
  
  # If conditions are not met, add the row to non_matched_rows
  if (!satisfied_conditions) {
    non_matched_rows <- rbind(non_matched_rows, row1)
  }
}
# # Write the matched_rows data frame to a BED file
# write.table(matched_rows, "matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
# Write the non_matched_rows data frame to a BED file
write.table(non_matched_rows, "non_matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

bedPeaksFile = "non_matched_rows.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
unmatched_column <- peakAnno_df[, 17]
matched_column <- matched_rows[, 5]

combined_column <- c(unmatched_column, matched_column)

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
unique_genes <- unique(combined_column)
gene_map <- getBM(
  filters = "mgi_symbol",          # Adjust if your gene names are not MGI symbols
  attributes = c("mgi_symbol", "ensembl_gene_id"),
  values = unique_genes,
  mart = mart
)
gene_map_vector <- setNames(gene_map$ensembl_gene_id, gene_map$mgi_symbol)
combined_column_converted <- gene_map_vector[combined_column]
combined_data <- data.frame(Value = combined_column_converted)
write.csv(combined_data, "FE_annotation.csv", row.names = FALSE)

#######################FG##################################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/FG")
Differential_FG = read.csv('Differential_FG.csv')
FG<- Differential_FG[,c(2,3,4,5)]
write.table(FG,'FG.bed',
            row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
file1 <- read.table("FG.bed", header = FALSE)
file2 <- read.csv("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/HiC/Cortex_HiC.csv", header = FALSE)

# Create an empty data frame to store matching rows
matched_rows <- data.frame()

# Create an empty data frame to store non-matching rows from file1
non_matched_rows <- data.frame()

# Iterate through each row in file1
for (i in 1:nrow(file1)) {
  # Get the current row from file1
  row1 <- file1[i, ]
  
  # Find matching rows in file2 based on column 1
  matching_rows_file2 <- file2[file2$V1 == row1$V1, ]
  
  # Initialize a flag to determine if the row1 satisfies the conditions
  satisfied_conditions <- FALSE
  
  # Iterate through matching rows in file2
  for (j in 1:nrow(matching_rows_file2)) {
    # Get the current matching row from file2
    row2 <- matching_rows_file2[j, ]
    
    # Check conditions
    if (row1$V2 <= row2$V3 && row1$V3 >= row2$V2) {
      # Set the flag to TRUE if conditions are met
      satisfied_conditions <- TRUE
      # Add the row from file1 to the data frame if conditions are met
      # Also add column 4 of file 2 to the new data frame
      matched_row <- c(row1, row2[4])  # Combine columns from both rows
      matched_rows <- rbind(matched_rows, matched_row)
      break  # Break out of the inner loop once a match is found
    }
  }
  
  # If conditions are not met, add the row to non_matched_rows
  if (!satisfied_conditions) {
    non_matched_rows <- rbind(non_matched_rows, row1)
  }
}
# # Write the matched_rows data frame to a BED file
# write.table(matched_rows, "matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
# Write the non_matched_rows data frame to a BED file
write.table(non_matched_rows, "non_matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

bedPeaksFile = "non_matched_rows.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
unmatched_column <- peakAnno_df[, 17]
matched_column <- matched_rows[, 5]

combined_column <- c(unmatched_column, matched_column)

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
unique_genes <- unique(combined_column)
gene_map <- getBM(
  filters = "mgi_symbol",          # Adjust if your gene names are not MGI symbols
  attributes = c("mgi_symbol", "ensembl_gene_id"),
  values = unique_genes,
  mart = mart
)
gene_map_vector <- setNames(gene_map$ensembl_gene_id, gene_map$mgi_symbol)
combined_column_converted <- gene_map_vector[combined_column]
combined_data <- data.frame(Value = combined_column_converted)
write.csv(combined_data, "FG_annotation.csv", row.names = FALSE)
