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
args = commandArgs(trailingOnly=TRUE)


deconv_results <- readRDS(file.path(here("temp_output", "cell_deconvolution", "deconv_results.rds")))
seurat_obj <- readRDS(here("temp_output", "preprocess", "filtered_data.rds"))

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
# find most abundant cell types
to_compare.cell_types <- cell_types_present[order(colSums(table(cell_type.df)), decreasing = TRUE)[1:2]]
# use only capture spots in those categories
# filtered.cell_type.df <- filter(cell_type.df, cell_type %in% to_compare.cell_types)

# group1 is the most abundant, group 2 is second most abundant
# group1 <- filter(cell_type.df, cell_type %in% to_compare.cell_types[1])%>%select(spot)
# group2 <- filter(cell_type.df, cell_type %in% to_compare.cell_types[2])%>%select(spot)

#'*step 3: find DEGs*
#'[IMPORTANT: positive values mean enrichment in GROUP1]
seurat_obj <- NormalizeData(seurat_obj)

#'*subset cells to run faster for testing
# Extract spatial coordinates
cd <- GetAssayData(seurat_obj, assay="Spatial", layer="counts") 
pos <- GetTissueCoordinates(seurat_obj)

# Randomly select 10 cells
random_cells <- sample(rownames(pos), 10)

# Subset Seurat object
seurat_obj <- subset(seurat_obj, cells = random_cells)

degs <- FindMarkers(seurat_obj, ident.1 = to_compare.cell_types[1], ident.2 = to_compare.cell_types[2], group.by = "cell_type", test="negbinom")

#TODO: volcano plot? 
# print(degs)

#'*step 4: save results*
savedir <- here("temp_output", "degs")
saveRDS(degs, file.path(savedir, "degs.rds"))

