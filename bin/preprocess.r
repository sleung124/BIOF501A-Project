#!/usr/bin/env Rscript

library(Seurat)
library(tidyverse)
library(here)

params.MITO_THRESH <- 20
params.path_to_sample <- here("data", "test")

# load in data, just from local directory for now
data <- Load10X_Spatial(params.path_to_sample)

# filter capture spots with high mitochondrial contamination
data <- PercentageFeatureSet(data, "^(MT|mt)-", col.name="percent.mito")
print(paste0("Found ", sum(data@meta.data$percent.mito > params.MITO_THRESH), " mitochondria-contaminated cells"))
data <- subset(data, subset = percent.mito <= params.MITO_THRESH)

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
mitoplot_path = here("temp_output", "preprocess")
saveRDS(data, file = file.path(mitoplot_path, "filtered_data.rds"))
ggsave(file.path(mitoplot_path, "mitoplot.jpg"), mitoplot, height = 8, width = 10)
