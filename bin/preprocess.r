#!/usr/bin/env Rscript

library(Seurat)
library(tidyverse)
library(here)

#'*commented out to test parameter initialization*
# MITO_THRESHOLD <- 20
# PATH_TO_SAMPLE <- here("data", "test")

# Use command line arguments. If user doesn't specify them, use 
# default parameters as outlined in nextflow.config

args = commandArgs(trailingOnly=TRUE)

MITO_THRESHOLD <- args[1]
PATH_TO_SAMPLE <- args[2]

# print(PATH_TO_SAMPLE)

# my_df <- data.frame(rep(10, 5), rep(20, 5))
# not_want_df <- my_df[, -1]

# write.csv(my_df, "dummy.csv")

# write.csv(not_want_df, "not_want.csv")

# load in data, just from local directory for now
if (grepl("temp_sample.rds", PATH_TO_SAMPLE)) { data <- readRDS(PATH_TO_SAMPLE)
} else {data <- Load10X_Spatial(PATH_TO_SAMPLE)}

# data <- readRDS(PATH_TO_SAMPLE)

# filter capture spots with high mitochondrial contamination
data <- PercentageFeatureSet(data, "^(MT|mt)-", col.name="percent.mito")
print(paste0("Found ", sum(data@meta.data$percent.mito > MITO_THRESHOLD), " mitochondria-contaminated cells"))
data <- subset(data, subset = percent.mito <= MITO_THRESHOLD)

# plot mitochondrial count distribution
mitoplot <- SpatialFeaturePlot(data, features="percent.mito",
                   alpha=0.7, pt.size.factor=1.2, 
                   crop=TRUE) + 
  ggtitle("") +
  scale_fill_gradientn(
    name = "Percentage",
    colors=pals::kovesi.rainbow(20)) + 
  theme_bw() +
  coord_fixed()  +
  theme(
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    legend.box.spacing = unit(2, "cm"),
    legend.key.size = unit(1, "cm"),
    legend.title = element_text(size=14, vjust=3.5),
    legend.text = element_text(size=14))

# save filtered data and mitoplot
# Create the folder if it doesn't exist
# folder_path <- here("temp_output")
# if (!dir.exists(folder_path)) {
#     dir.create(folder_path)
#     folders <- c("preprocess", "cell_deconvolution", "pathways", "degs")
#     for (i in 1:length(folders)) {
#         dir.create(here(folder_path, folders[i]))
#     }
# } 

# mitoplot_path = here("temp_output", "preprocess")
# saveRDS(data, file = file.path(mitoplot_path, "filtered_data.rds"))
saveRDS(data, file = "filtered_data.rds")
ggsave("mitoplot.jpg", mitoplot, height = 8, width = 10)
