############preparation##############
##This is just carry-over, maybe delete###
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
####################################BA1_D1C_enhancer###############################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/HiC")
# Read the two input files
##If there are diff peaks on chrY, add in row number for chry
file1 <- read.table("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/overlap/BA1_D1C/overlapped_BA1_D1C_enhancer.bed", header = FALSE)
file2 <- read.csv("Cortex_HiC.csv", header = FALSE)
# file1 <- file11[-2853,]
#file1 <- file111[-1301,]

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
write.table(non_matched_rows, "overlapped_BA1_D1C_non_matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

bedPeaksFile = "overlapped_BA1_D1C_non_matched_rows.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
anno = grepl('Promoter', peakAnno_df$annotation,fixed = TRUE)
new_peakAnno_df = peakAnno_df[!anno,]
unmatched_column <- new_peakAnno_df[, 17]
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
write.csv(combined_data, "overlapped_BA1_D1C_enhancerHiC_annotation.csv", row.names = FALSE)

####################################CA1_BD1_enhancer###############################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/HiC")
# Read the two input files
##If there are diff peaks on chrY, add in row number for chry
file1 <- read.table("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/overlap/CA1_BD1/overlapped_CA1_BD1_enhancer.bed", header = FALSE)
file2 <- read.csv("Cortex_HiC.csv", header = FALSE)
# file1 <- file11[-2853,]
#file1 <- file111[-1301,]

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
write.table(non_matched_rows, "overlapped_CA1_BD1_non_matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

bedPeaksFile = "overlapped_CA1_BD1_non_matched_rows.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
anno = grepl('Promoter', peakAnno_df$annotation,fixed = TRUE)
new_peakAnno_df = peakAnno_df[!anno,]
unmatched_column <- new_peakAnno_df[, 17]
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
write.csv(combined_data, "overlapped_CA1_BD1_enhancerHiC_annotation.csv", row.names = FALSE)

####################################BA1_UNIQUE_enhancer###############################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/HiC")
# Read the two input files
##If there are diff peaks on chrY, add in row number for chry
file1 <- read.table("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/overlap/subtract_BA1_enhancer.bed", header = FALSE)
file2 <- read.csv("Cortex_HiC.csv", header = FALSE)
# file1 <- file11[-2853,]
#file1 <- file111[-1301,]

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


write.table(non_matched_rows, "subtract_BA1_enhancer_non_matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

bedPeaksFile = "subtract_BA1_enhancer_non_matched_rows.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
anno = grepl('Promoter', peakAnno_df$annotation,fixed = TRUE)
new_peakAnno_df = peakAnno_df[!anno,]
unmatched_column <- new_peakAnno_df[, 17]
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

write.csv(combined_data, "subtract_BA1_enhancerHiC_annotation.csv", row.names = FALSE)



####################################D1C_UNIQUE_enhancer###############################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/HiC")
# Read the two input files
##If there are diff peaks on chrY, add in row number for chry
file1 <- read.table("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/overlap/subtract_D1C_enhancer.bed", header = FALSE)
file2 <- read.csv("Cortex_HiC.csv", header = FALSE)
# file1 <- file11[-2853,]
#file1 <- file111[-1301,]

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


write.table(non_matched_rows, "subtract_D1C_enhancer_non_matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

bedPeaksFile = "subtract_D1C_enhancer_non_matched_rows.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
anno = grepl('Promoter', peakAnno_df$annotation,fixed = TRUE)
new_peakAnno_df = peakAnno_df[!anno,]
unmatched_column <- new_peakAnno_df[, 17]
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

write.csv(combined_data, "subtract_D1C_enhancerHiC_annotation.csv", row.names = FALSE)

####################################CA1_UNIQUE_enhancer###############################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/HiC")
# Read the two input files
##If there are diff peaks on chrY, add in row number for chry
file1 <- read.table("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/overlap/subtract_CA1_enhancer.bed", header = FALSE)
file2 <- read.csv("Cortex_HiC.csv", header = FALSE)
# file1 <- file11[-2853,]
#file1 <- file111[-1301,]

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

write.table(non_matched_rows, "subtract_CA1_enhancer_non_matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

bedPeaksFile = "subtract_CA1_enhancer_non_matched_rows.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
anno = grepl('Promoter', peakAnno_df$annotation,fixed = TRUE)
new_peakAnno_df = peakAnno_df[!anno,]
unmatched_column <- new_peakAnno_df[, 17]
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
write.csv(combined_data, "subtract_CA1_enhancerHiC_annotation.csv", row.names = FALSE)

####################################BD1_UNIQUE_enhancer###############################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/HiC")
# Read the two input files
##If there are diff peaks on chrY, add in row number for chry
file1 <- read.table("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/overlap/subtract_BD1_enhancer.bed", header = FALSE)
file2 <- read.csv("Cortex_HiC.csv", header = FALSE)
# file1 <- file11[-2853,]
#file1 <- file111[-1301,]

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

write.table(non_matched_rows, "subtract_BD1_enhancer_non_matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

bedPeaksFile = "subtract_BD1_enhancer_non_matched_rows.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
anno = grepl('Promoter', peakAnno_df$annotation,fixed = TRUE)
new_peakAnno_df = peakAnno_df[!anno,]
unmatched_column <- new_peakAnno_df[, 17]
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
write.csv(combined_data, "subtract_BD1_enhancerHiC_annotation.csv", row.names = FALSE)

####################################BA1_FE_enhancer###############################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/HiC")
# Read the two input files
##If there are diff peaks on chrY, add in row number for chry
file1 <- read.table("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/overlap/subtract_BA1_FE_enhancer.bed", header = FALSE)
file2 <- read.csv("Cortex_HiC.csv", header = FALSE)
# file1 <- file11[-2853,]
#file1 <- file111[-1301,]

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
write.table(non_matched_rows, "overlapped_BA1_D1C_non_matched_rows.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

bedPeaksFile = "overlapped_BA1_D1C_non_matched_rows.bed"
peak <- readPeakFile(bedPeaksFile)
keepChr=!grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak,
                         tssRegion = c(-2000, 2000),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")
peakAnno_df = as.data.frame(peakAnno)
anno = grepl('Promoter', peakAnno_df$annotation,fixed = TRUE)
new_peakAnno_df = peakAnno_df[!anno,]
unmatched_column <- new_peakAnno_df[, 17]
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
write.csv(combined_data, "subtract_BA1_FE_enhancer_annotation.csv", row.names = FALSE)
