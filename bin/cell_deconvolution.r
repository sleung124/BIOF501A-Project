#!/usr/bin/env Rscript

# Rscript for reference-free cell deconvolution

library(STdeconvolve)
library(Seurat)
library(tidyverse)
library(here)

# Parameters from user or config file
args <- commandArgs(trailingOnly=TRUE)

MAX_LDA_K <- as.integer(args[1])
RADIUS <- as.double(args[2])
data <- readRDS(args[3])

# adapting start code from STdeconvolve docs, link below:
# https://github.com/JEFworks-Lab/STdeconvolve
cd <- GetAssayData(data, assay="Spatial", layer="counts") 
pos <- GetTissueCoordinates(data)
colnames(pos) <- c("y", "x")
pos$y <- -1*pos$y

# remove pixels with low amount of genes
counts <- cleanCounts(cd, min.lib.size=100, min.reads=10)
# feature select for genes
corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow=0.05)
## choose optimal number of cell-types
# this takes a really really long time; default use just 2
ldas <- fitLDA(t(as.matrix(corpus)), Ks = seq(2, MAX_LDA_K, by = 1))
## get best model results
optLDA <- optimalModel(models = ldas, opt = "min")
## extract deconvolved cell-type proportions (theta) and transcriptional profiles (beta)
results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
deconProp <- results$theta
deconGexp <- results$beta

## visualize deconvolved cell-type proportions
plt <- vizAllTopics(theta = deconProp,
                    pos = pos,
                    r = RADIUS,
                    lwd = 0,
                    showLegend = TRUE,
                    plotTitle = NA) + 
  ggplot2::guides(fill=ggplot2::guide_legend(ncol=2)) +
  ## outer border
  ggplot2::geom_rect(data = data.frame(pos),
                     ggplot2::aes(xmin = min(x)-90, xmax = max(x)+90,
                                  ymin = min(y)-90, ymax = max(y)+90),
                     fill = NA, color = "black", linetype = "solid", linewidth = 0.5) +
  ggplot2::theme(
    plot.background = ggplot2::element_blank()
  ) +
  ## remove the pixel "groups", which is the color aesthetic for the pixel borders
  ggplot2::guides(colour = "none")

# save generated plot and deconvolution results
saveRDS(results, file = "deconv_results.rds")
ggsave("deconvolution.jpg", plt, height = 8, width = 10)
