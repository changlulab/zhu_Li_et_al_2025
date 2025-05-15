############ Preparation ############
library(ggplot2)
library(dplyr)
library(stats)
library(gplots)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

############ Creating your new count file containing just diff peaks ############
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/count_files")
dataframe_A <- read.csv("enhancer_count_abcd.csv")

# Read BED files
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer/old files")
bed_files <- lapply(c("Diff_BA1_peaks.bed", "Diff_D1C_peaks.bed", "Diff_CA1_peaks.bed", 
                      "Diff_BD1_peaks.bed", "Diff_D1A1_peaks.bed"), read.table,
                    header = FALSE, stringsAsFactors = FALSE)

dataframe_B <- do.call(rbind, bed_files)
values_to_match <- unique(dataframe_B[, 2])
matching_rows <- dataframe_A[dataframe_A[, 2] %in% values_to_match, ]
filtered_columns <- matching_rows %>% select(1:15, 33:68)

write.csv(filtered_columns, "all_counts.csv", row.names = FALSE)

# Normalize using L2 norm
new_dataframe <- filtered_columns %>% select(5:ncol(filtered_columns))
normalize_l2 <- function(x) x / sqrt(sum(x^2))
normalized_rows <- t(apply(new_dataframe, 1, normalize_l2))
normalized_dataframe <- as.data.frame(normalized_rows)

# Extract groups
GroupA_mean <- rowMeans(normalized_dataframe %>% select(matches("Group.A")), na.rm = TRUE)
GroupB_mean <- rowMeans(normalized_dataframe %>% select(matches("Group.B")), na.rm = TRUE)
GroupC_mean <- rowMeans(normalized_dataframe %>% select(matches("Group.C")), na.rm = TRUE)
GroupD_mean <- rowMeans(normalized_dataframe %>% select(matches("Group.D")), na.rm = TRUE)
row_means <- rowMeans(cbind(GroupA_mean, GroupB_mean, GroupC_mean, GroupD_mean), na.rm = TRUE)

input_1 <- data.frame(
  START = filtered_columns[, 2],
  GroupA = GroupA_mean,
  GroupB = GroupB_mean,
  GroupC = GroupC_mean,
  GroupD = GroupD_mean,
  Average = row_means
)

input_2 <- data.frame(
  CHR = filtered_columns[, 1],
  START = filtered_columns[, 2],
  END = filtered_columns[, 3],
  GroupA = GroupA_mean,
  GroupB = GroupB_mean,
  GroupC = GroupC_mean,
  GroupD = GroupD_mean,
  Average = row_means
)

write.csv(input_1, "ABCD_input_norm.csv", row.names = FALSE)
write.csv(input_2, "ABCD_input_norm_2.csv", row.names = FALSE)

############ K-means for 4 groups ############
my_data1 <- read.csv("ABCD_input_norm.csv")
rownames(my_data1) <- make.unique(as.character(my_data1$START))  # Ensure unique rownames

my_data <- my_data1[, 2:5]

set.seed(123)
k <- 4
kmeans_result <- kmeans(my_data, centers = k, nstart = 23)
cluster_labels <- kmeans_result$cluster

# Save genes in each cluster
for (cluster in 1:k) {
  genes_in_cluster <- rownames(my_data1)[cluster_labels == cluster]
  write.csv(genes_in_cluster, file = paste0("cluster_", cluster, "_genes_ABCD.csv"), row.names = FALSE)
}

# Add cluster info and prepare for heatmaps
my_data2 <- read.csv("ABCD_input_norm_2.csv")
my_data2$Cluster <- cluster_labels
my_data2 <- my_data2[order(my_data2$Cluster), ]
mydata3 <- as.matrix(my_data2[, 4:7])

# Generate labels
unique_labels <- unique(my_data2$Cluster)
labels_row <- rep("", nrow(mydata3))
for (i in 1:length(unique_labels)) {
  indices <- which(my_data2$Cluster == unique_labels[i])
  labels_row[indices[1]] <- paste("Cluster", as.character(unique_labels[i]))
}

############ pheatmap: No clustering ############
pdf(file = "kmeans_enhancer_average_color_change.pdf", width = 4.5, height = 7)
pheatmap(mydata3,
         scale = "none",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         labels_row = labels_row,
         border_color = NA,
         fontsize_col = 12,
         fontsize_row = 10,
         fontsize = 8,
         color = colorRampPalette(c("white", "#D6F4FF", "#0c84c6"))(100),
         angle_col = "0")
dev.off()

############ ComplexHeatmap: Bottom legend with black border and black ticks ############
col_fun <- colorRamp2(seq(min(mydata3), max(mydata3), length.out = 100),
                      colorRampPalette(c("white", "#D6F4FF", "#0c84c6"))(100))

ht <- Heatmap(mydata3,
              name = "My Heatmap",
              col = col_fun,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_row_names = TRUE,
              show_column_names = TRUE,
              row_names_gp = gpar(fontsize = 10),
              column_names_gp = gpar(fontsize = 12),
              row_title = "Clusters",
              column_title = "Groups",
              row_labels = labels_row,
              heatmap_legend_param = list(
                legend_direction = "horizontal",
                title = "Scale",
                title_position = "topcenter",
                legend_width = unit(6, "cm"),
                border = "black",               # Black border around legend
                ticks_gp = gpar(col = "black")  # Black tick marks
              ))

pdf(file = "kmeans_enhancer_average_color_change_scale.pdf", width = 5, height = 4)
draw(ht, heatmap_legend_side = "bottom")
dev.off()

############ pheatmap: Column clustering ON ############
pdf(file = "kmeans_enhancer_column_cluster_color_change_average.pdf", width = 5, height = 4)
pheatmap(mydata3,
         scale = "none",
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         labels_row = labels_row,
         border_color = NA,
         fontsize_row = 8,
         fontsize = 8,
         color = colorRampPalette(c("white", "#D6F4FF", "#0c84c6"))(100))
dev.off()


############Kmeans_clustering_gene_annotation##########################################################
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
library(readr) # Ensure to include readr for read_csv function
options(stringsAsFactors = F)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/ChIP/result/differential peaks/enhancer")

# Function to extract and save data for each cluster
extract_and_save <- function(cluster_number) {
  # Construct file paths
  cluster_file <- paste0("cluster_", cluster_number, "_genes_ABCD.csv")
  output_file <- paste0("cluster_", cluster_number, ".bed")
  
  # Read row numbers to extract
  rows_to_extract <- read_csv(cluster_file, col_names = FALSE, skip = 1)
  rows_to_extract <- unlist(rows_to_extract)
  
  # Read all_counts.csv
  all_counts <- read_csv("all_counts.csv")
  
  # Extract rows based on row numbers
  extracted_data <- all_counts[rows_to_extract, , drop = FALSE]
  
  # Save first three columns of extracted data to a BED file
  write.table(extracted_data[, 1:3], output_file, col.names = FALSE,
              row.names = FALSE, sep = "\t", quote = FALSE)
}

# Loop through clusters 1 to 8 to extract and save BED files
for (cluster_number in 1:4) {
  extract_and_save(cluster_number)
}

# Function to process cluster data for annotation
process_cluster_annotation <- function(cluster_number) {
  # Read the BED file for the cluster
  bed_file <- paste0("cluster_", cluster_number, ".bed")
  file1 <- read.table(bed_file, header = FALSE)
  
  # Read the HiC data
  file2 <- read.csv("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/MIA data/HiC/Cortex_HiC.csv", header = FALSE)
  
  # Create empty data frames for matching and non-matching rows
  matched_rows <- data.frame()
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
        # Also add column 4 of file2 to the new data frame
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
  
  # Write the non_matched_rows data frame to a BED file
  non_match_file <- paste0("cluster_", cluster_number, "_non_match.bed")
  write.table(non_matched_rows, non_match_file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  # Read and annotate the non-matching peaks
  bedPeaksFile <- non_match_file
  peak <- readPeakFile(bedPeaksFile)
  keepChr <- !grepl('_', seqlevels(peak))
  seqlevels(peak, pruning.mode = "coarse") <- seqlevels(peak)[keepChr]
  peakAnno <- annotatePeak(peak,
                           tssRegion = c(-2000, 2000),
                           TxDb = txdb,
                           annoDb = "org.Mm.eg.db")
  peakAnno_df <- as.data.frame(peakAnno)
  
  # Extract unmatched and matched columns
  unmatched_column <- peakAnno_df[, 16]
  matched_column <- matched_rows[, 4]
  
  combined_column <- c(unmatched_column, matched_column)
  
  # Retrieve gene mappings using biomaRt
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
  
  # Write the combined data to a CSV file
  combined_data <- data.frame(Value = combined_column_converted)
  output_csv <- paste0("cluster_", cluster_number, "_annotation.csv")
  write.csv(combined_data, output_csv, row.names = FALSE)
}

# Loop through clusters 1 to 8 for annotation
for (cluster_number in 1:4) {
  process_cluster_annotation(cluster_number)
}