#!/usr/bin/env Rscript


# load libaries and default params
library(enrichR)
library(tidyverse)

#'*Nextflow Params*
listed_dbs <- listEnrichrDbs()
# user can provide a list of libraries to look through 
DBS <- unique(listed_dbs$libraryName)[grep("Mouse",unique(listed_dbs$libraryName))]
PVAL_THRESH <- 0.05
# params for plotEnrich function
SHOW_TERMS <- 20
NUMCHAR <- 40
ORDER_BY <- "P.value"

degs <- readRDS(here("temp_output", "degs", "degs.rds")) %>%
  filter(p_val_adj < PVAL_THRESH & p_val_adj > 0) %>%
  rownames_to_column("genes") %>%
  select(genes) %>%
  pull()

enriched <- enrichr(degs, params.DBS)

#'*Save plots of enriched pathways*
savedir <- here("temp_output", "pathways")
for (i in 1:length(enriched)) {
  plt_name <- file.path(savedir, paste0(names(enriched[i]), "pathways.jpg"))
  plt <- plotEnrich(enriched[[i]], showTerms = SHOW_TERMS, numChar = NUMCHAR, orderBy = ORDER_BY)
  ggsave(plt_name, plt)
}

