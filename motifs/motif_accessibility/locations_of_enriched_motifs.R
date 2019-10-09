library(data.table)
library(purrr)


io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$motif_bed <- paste0(io$data_dir, "/features/motifs/260419/200bp")
io$sig_motifs <- "/bi/home/clarks/gastrulation/plots/motifs/diffexp_diffacc/background_all_k27ac/260419/up"

io$out <- paste0(io$data_dir, "/features/lineage_enriched_motif_positions/")

opts<-list()
opts$dif_acc <- c(#"Ectoderm_sites_TFs_EPIvsEct",
                  "E7.5Mes_vs_E7.5EctEnd_H3K27ac_distal_E7.5_Mes_intersect12",
                  "E7.5End_vs_E7.5EctMes_H3K27ac_distal_E7.5_End_intersect12",
                  "E7.5Ect_vs_E7.5EndMes_H3K27ac_distal_E7.5_Ect_intersect12")



sig_motifs <- dir(io$sig_motifs, pattern = ".tsv", full = TRUE) %>%
  set_names(basename(.) %>% gsub(".tsv", "", .)) %>%
  .[names(.) %in% opts$dif_acc] %>%
  map(fread) %>%
  map(~.[sig_both == TRUE, .(tf, motif, comparison_rna = comparison, anno)]) %>%
  map2(names(.), ~.x[, comparison_acc := .y]) %>%
  rbindlist() %>%
  split(by = c("comparison_rna", "comparison_acc"))


motif_positions <- dir(io$motif_bed, full = TRUE) %>%
  set_names(basename(.) %>% gsub(".bed", "", .)) %>%
  map(fread)

bed_files <- map(sig_motifs, ~{
  an <- .[, unique(anno)] %>%
    .[1]
  tfs <- .[, unique(motif)]
  motif_positions[names(motif_positions) %like% an] %>%
    rbindlist() %>%
    .[tf %in% tfs] %>%
    split(by = "tf")
})




