#!/usr/bin/env Rscript


# load libaries and default params
library(enrichR)
library(tidyverse)
library(here)

#'* declare params from config file / user given arguments
args = commandArgs(trailingOnly=TRUE)

PVAL_THRESH <- as.double(args[1])
SHOW_TERMS <- as.integer(args[2])
NUMCHAR <- as.integer(args[3])
ORDER_BY <- args[4]
SPECIES <- paste0(toupper(substr(args[5], 1, 1)), substr(args[5], 2, nchar(args[5])))
LOADED_DEGS <- readRDS(args[6])

# TODO: print out args to debug; think SPECIES is not parsing properly
print(paste0("Argument for processed Species: ", SPECIES))

#'*Nextflow Params*
listed_dbs <- listEnrichrDbs()
# user can provide a list of libraries to look through 

DBS <- unique(listed_dbs$libraryName)[grep(SPECIES,unique(listed_dbs$libraryName))]

degs <- LOADED_DEGS %>%
  filter(p_val_adj < PVAL_THRESH & p_val_adj > 0) %>%
  rownames_to_column("genes") %>%
  select(genes) %>%
  pull()

enriched <- enrichr(degs, DBS)

#'*Save plots of enriched pathways*
for (i in 1:length(enriched)) {
  print(paste0("On output ", i, " out of ", length(enriched)))
  plt_name <- file.path(paste0(names(enriched[i]), "_pathways.jpg"))
  plt <- plotEnrich(enriched[[i]], showTerms = SHOW_TERMS, numChar = NUMCHAR, orderBy = ORDER_BY)
  ggsave(plt_name, plt)
}

