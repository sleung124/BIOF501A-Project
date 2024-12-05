#!/usr/bin/env Rscript

# use seurat FindMarkers for DEGs

#'*general steps for now:* 
# 1) load in results from cell deconv
# 2) group each capture spot by most abundant cell type
# 3) find most abundant cell types and perform deg analysis

# step 0: load in libraries needed; declare params
library(here)
library(tidyverse)
library(Seurat)

#'*step 1: load in data*
args <- commandArgs(trailingOnly=TRUE)

quick_sample <- as.integer(args[1])
seed <- as.integer(args[2])
PVAL_THRESH <- as.double(args[3])
deconv_results <- readRDS(args[4])
seurat_obj <- readRDS(args[5])

if (seed > 0) {
  set.seed(seed)
}
#'*step 2: categorize capture spots by cells*
# data$theta is a matrix. Rows are capture spots, columns are cell type. Intersection is cell proportion
cell.proportions <- deconv_results$theta
number_capture_spots <- dim(cell.proportions)[1]
cell_types_present <- colnames(cell.proportions)
primary_type <- rep(NA, number_capture_spots)

for (i in 1:number_capture_spots) {
  primary_type[i] <- cell_types_present[which.max(cell.proportions[i,])]
}

# assign predicted majority cell type to seurat object
seurat_obj@meta.data$cell_type <- primary_type

cell_type.df <- data.frame(spot = rownames(cell.proportions), cell_type = primary_type)

# find most abundant cell types, compare the two most abundant cell types
# use only capture spots from these regions
to_compare.cell_types <- cell_types_present[order(colSums(table(cell_type.df)), decreasing = TRUE)[1:2]]

#'*step 3: find DEGs*
#'[IMPORTANT: positive values mean enrichment in GROUP1, which is the most abundant cell type group]
seurat_obj <- NormalizeData(seurat_obj)

#'*subset cells to run faster for testing
# samples for barcode names of capture spots
# ONLY USED TO SPEED UP TEST EXAMPLE, SHOULD NOT USE FOR REAL DATA ANALYSIS
if (quick_sample > 0) {
  print(paste0("Sampling ", quick_sample, " random capture spots"))

  cd <- GetAssayData(seurat_obj, assay="Spatial", layer="counts") 
  pos <- GetTissueCoordinates(seurat_obj)

  # Randomly select 10 cells
  random_cells <- sample(rownames(pos), quick_sample)

  # Subset Seurat object
  seurat_obj <- subset(seurat_obj, cells = random_cells)
}

# find DEGs
degs <- FindMarkers(seurat_obj, ident.1 = to_compare.cell_types[1], ident.2 = to_compare.cell_types[2], group.by = "cell_type", test="negbinom") 

# Generate volcano plot
volcano_data <- degs %>%
  mutate(
    Groups = ifelse(p_val_adj < PVAL_THRESH & avg_log2FC > 0, "Most Abundant Cell Type",
                   ifelse(p_val_adj < PVAL_THRESH & avg_log2FC < 0, "2nd Most Abundant Cell Type", 
                          "Non-Significant"))
  )

volcano_plot <- volcano_data %>%
  ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj), color = Groups)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("Most Abundant Cell Type" = "red", 
                                "2nd Most Abundant Cell Type" = "blue",
                                "Not Significant" = "grey")) +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-Value",
  ) +
  theme_minimal() +
  ggtitle("Volcano Plot - Cell Type Abundance Comparisons")

#'*step 4: save results*
sig_degs <- degs %>% filter(p_val_adj < PVAL_THRESH & p_val_adj > 0)
write.csv(degs, "degs.csv")
ggsave("volcano_plot.jpg", volcano_plot, height = 8, width = 10)
saveRDS(degs, "degs.rds")
