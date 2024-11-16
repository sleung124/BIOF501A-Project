#!/usr/bin/env Rscript

library(Seurat)
library(here)

# running into issue listed here: https://github.com/satijalab/seurat/issues/9061

to_save <- c("anterior1", "anterior2", "posterior1", "posterior2")
for (i in 1:length(to_save)) {
  file_name <- paste0(to_save[i], ".rds")
  assign(to_save[i], readRDS(file = here("temp", "data", file_name)))
}

plot1 <- VlnPlot(anterior1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot1.1 <- SpatialFeaturePlot(anterior1, features = "nCount_Spatial")+ NoLegend()

plot2 <- VlnPlot(anterior2, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2.1 <- SpatialFeaturePlot(anterior2, features = "nCount_Spatial")+ NoLegend()

plot3 <- VlnPlot(posterior1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot3.1 <- SpatialFeaturePlot(posterior1, features = "nCount_Spatial")+ NoLegend()

plot4 <- VlnPlot(posterior2, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot4.1 <- SpatialFeaturePlot(posterior2, features = "nCount_Spatial")+ NoLegend()

wrap_plots(plot1.1, plot2.1, plot3.1, plot4.1)


