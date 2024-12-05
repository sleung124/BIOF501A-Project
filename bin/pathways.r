#!/usr/bin/env Rscript

# Rscript dedicated to finding the enriched pathways present
# Uses output of find_degs.r

# load libaries and default params
library(enrichR)
library(tidyverse)
library(here)

#'* declare params from config file / user given arguments
args <- commandArgs(trailingOnly=TRUE)

# PVAL_THRESH <- as.double(args[1])
SHOW_TERMS <- as.integer(args[1])
NUMCHAR <- as.integer(args[2])
ORDER_BY <- args[3]
SPECIES <- paste0(toupper(substr(args[4], 1, 1)), substr(args[4], 2, nchar(args[4])))
LOADED_DEGS <- readRDS(args[5])

# List of possible databases 
listed_dbs <- listEnrichrDbs()

# Filter for Mouse or Human databases
if (SPECIES == "Mouse") {
  cond <- grep(SPECIES,unique(listed_dbs$libraryName))
} else {
  cond <- !grep(SPECIES,unique(listed_dbs$libraryName))
}
DBS <- unique(listed_dbs$libraryName)[cond]

# Grab only gene names
degs <- LOADED_DEGS  %>%
  rownames_to_column("genes") %>%
  select(genes) %>%
  pull()

# Query enrichr
enriched <- enrichr(degs, DBS)

#'*Save plots of enriched pathways*
for (i in 1:length(enriched)) {
  print(paste0("On output ", i, " out of ", length(enriched)))
  plt_name <- file.path(paste0(names(enriched[i]), "_pathways.jpg"))
  plt <- plotEnrich(enriched[[i]], showTerms = SHOW_TERMS, numChar = NUMCHAR, orderBy = ORDER_BY)
  ggsave(plt_name, plt)
}

