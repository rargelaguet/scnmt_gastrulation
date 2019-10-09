
library(data.table)
library(purrr)
library(ggplot2)
library(cowplot)


io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/EB_data_new_May2019"
io$in_dir <- "/bi/scratch/Stephen_Clark/EB_data_new_May2019/acc/differential/feature_level/motif_enrichment/background_all_k27ac"



motifs <- dir(io$in_dir, pattern = "ame.txt", recursive = TRUE, full = TRUE) %>%
  set_names(basename(dirname(.))) %>%
  map(fread) %>%
  map2(names(.), ~.x[, dif_met := .y]) %>%
  rbindlist() %>%
  .[, .(dif_met, tf = V9, p = V13, q = as.numeric(gsub(")", "", V16)) )]
