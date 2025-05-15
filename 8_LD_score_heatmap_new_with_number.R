rm(list = ls())
options(stringsAsFactors = FALSE)

# Load libraries
library(ComplexHeatmap)
library(circlize)
library(grid)

# Load data
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures/")
enhancer <- read.table("enhancer_ldsc_1.txt", header = TRUE, row.names = 1, sep = "\t")

# Order rows by average –log(P) value (descending)
row_means <- rowMeans(enhancer, na.rm = TRUE)
enhancer <- enhancer[order(-row_means), ]

# Split into panels
heatmap1 <- enhancer[, 1:4]
heatmap2 <- enhancer[, 5:7]

# Custom color palette using smooth interpolation
custom_breaks <- c(0, 1, 3, 5, 7, 10, 20)
base_colors <- c("white", "#e0f3ff", "dodgerblue1", "dodgerblue3", "blue3", "blue4")
smooth_colors <- colorRampPalette(base_colors)(length(custom_breaks))
color_blues <- colorRamp2(custom_breaks, smooth_colors)

# Shared legend settings
legend_param <- list(
  title = "-log(P)",
  direction = "horizontal",
  title_position = "topcenter",
  legend_width = unit(3, "cm"),
  border = "black",
  at = custom_breaks
)

# Cell drawing function with adaptive text color
make_clean_cell_fun <- function(mat, col_fun) {
  function(j, i, x, y, width, height, fill) {
    val <- mat[i, j]
    grid.rect(
      x = x, y = y,
      width = width, height = height,
      gp = gpar(fill = col_fun(val), col = NA)
    )
    if (!is.na(val)) {
      text_col <- if (val > 8) "white" else "black"
      grid.text(
        label = sprintf("%.1f", val),
        x = x, y = y,
        gp = gpar(fontsize = 8, col = text_col)
      )
    }
  }
}

# Panel 1
ht1 <- Heatmap(
  heatmap1,
  name = "-log(P)",
  col = color_blues,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10),
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 10),
  column_names_side = "top",
  column_names_rot = 0,
  column_names_centered = TRUE,
  rect_gp = gpar(col = NA),
  cell_fun = make_clean_cell_fun(heatmap1, color_blues),
  border = TRUE,
  heatmap_legend_param = legend_param
)

# Panel 2
ht2 <- Heatmap(
  heatmap2,
  name = "-log(P)",
  col = color_blues,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 10),
  column_names_side = "top",
  column_names_rot = 0,
  column_names_centered = TRUE,
  rect_gp = gpar(col = NA),
  cell_fun = make_clean_cell_fun(heatmap2, color_blues),
  border = TRUE
)

# Save to PDF
pdf("LD_score_enhancer_FINAL_with_values_adaptive_text.pdf", width = 6, height = 3.5)
draw(ht1 + ht2, heatmap_legend_side = "bottom", merge_legend = TRUE, gap = unit(0, "mm"))
dev.off()


##############promoter###################
rm(list = ls())
options(stringsAsFactors = FALSE)

# Load libraries
library(ComplexHeatmap)
library(circlize)
library(grid)

# Load data
setwd("/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures/")
promoter <- read.table("promoter_ldsc_1.txt", header = TRUE, row.names = 1, sep = "\t")

# Order rows by average –log(P) value (descending)
row_means <- rowMeans(promoter, na.rm = TRUE)
promoter <- promoter[order(-row_means), ]

# Split into panels
panel1 <- promoter[, 1:4]
panel2 <- promoter[, 5:ncol(promoter)]

# Compute value range
value_range <- range(promoter, na.rm = TRUE)
min_val <- floor(value_range[1])
max_val <- ceiling(value_range[2])

# Generate auto-scaled color gradient
auto_blues <- colorRamp2(
  seq(min_val, max_val, length.out = 7),
  colorRampPalette(c("white", "dodgerblue1", "blue4"))(7)
)

# Legend
legend_param <- list(
  title = "-log(P)",
  direction = "horizontal",
  title_position = "topcenter",
  legend_width = unit(3, "cm"),
  border = "black",
  at = pretty(seq(min_val, max_val, length.out = 7))
)

# Cell drawing function with text
make_clean_cell_fun <- function(mat, col_fun) {
  function(j, i, x, y, width, height, fill) {
    val <- mat[i, j]
    grid.rect(
      x = x, y = y,
      width = width, height = height,
      gp = gpar(fill = col_fun(val), col = NA)
    )
    if (!is.na(val)) {
      grid.text(
        label = sprintf("%.1f", val),
        x = x, y = y,
        gp = gpar(fontsize = 8, col = "black")
      )
    }
  }
}

# Panel 1
ht1 <- Heatmap(
  panel1,
  name = "-log(P)",
  col = auto_blues,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10),
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 10),
  column_names_side = "top",
  column_names_rot = 0,
  column_names_centered = TRUE,
  rect_gp = gpar(col = NA),
  cell_fun = make_clean_cell_fun(panel1, auto_blues),
  border = TRUE,
  heatmap_legend_param = legend_param
)

# Panel 2
ht2 <- Heatmap(
  panel2,
  name = "-log(P)",
  col = auto_blues,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 10),
  column_names_side = "top",
  column_names_rot = 0,
  column_names_centered = TRUE,
  rect_gp = gpar(col = NA),
  cell_fun = make_clean_cell_fun(panel2, auto_blues),
  border = TRUE
)

# Save to PDF
pdf("LD_score_promoter_FINAL_with_values.pdf", width = 6, height = 3.5)
draw(ht1 + ht2, heatmap_legend_side = "bottom", merge_legend = TRUE, gap = unit(0, "mm"))
dev.off()
