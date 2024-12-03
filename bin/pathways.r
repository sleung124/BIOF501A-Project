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
wanted_dbs <- args[5]
loaded_degs <- readRDS(args[6])

#'*Nextflow Params*
listed_dbs <- listEnrichrDbs()
# user can provide a list of libraries to look through 
if (wanted_dbs == "default_mouse") {
  DBS <- unique(listed_dbs$libraryName)[grep("Mouse",unique(listed_dbs$libraryName))]
} else {
  DBS <- wanted_dbs
}

degs <- loaded_degs %>%
  filter(p_val_adj < PVAL_THRESH & p_val_adj > 0) %>%
  rownames_to_column("genes") %>%
  select(genes) %>%
  pull()

enriched <- enrichr(degs, DBS)

#'*Save plots of enriched pathways*
# savedir <- here("temp_output", "pathways")
for (i in 1:length(enriched)) {
  plt_name <- file.path(paste0(names(enriched[i]), "pathways.jpg"))
  plt <- plotEnrich(enriched[[i]], showTerms = SHOW_TERMS, numChar = NUMCHAR, orderBy = ORDER_BY)
  ggsave(plt_name, plt)
}

