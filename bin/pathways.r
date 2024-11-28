#!/usr/bin/env Rscript


# load libaries and default params
library(enrichR)
library(tidyverse)

#'*Nextflow Params*
dbs <- listEnrichrDbs()
# user can provide a list of libraries to look through 
params.DBS <- unique(dbs$libraryName)[grep("Mouse",unique(dbs$libraryName))]
params.PVAL_THRESH <- 0.05
# params for plotEnrich function
params.SHOW_TERMS <- 20
params.NUMCHAR <- 40
params.ORDER_BY <- "P.value"

degs <- readRDS(here("temp_output", "degs", "degs.rds")) %>%
  filter(p_val_adj < params.PVAL_THRESH & p_val_adj > 0) %>%
  rownames_to_column("genes") %>%
  select(genes) %>%
  pull()

enriched <- enrichr(degs, params.DBS)

#'*Save plots of enriched pathways*
savedir <- here("temp_output", "pathways")
for (i in 1:length(enriched)) {
  plt_name <- file.path(savedir, paste0(names(enriched[i]), "pathways.jpg"))
  plt <- plotEnrich(enriched[[i]], showTerms = params.SHOW_TERMS, numChar = params.NUMCHAR, orderBy = params.ORDER_BY)
  ggsave(plt_name, plt)
}

