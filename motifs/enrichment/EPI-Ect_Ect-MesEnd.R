library(data.table)
library(purrr)
library(furrr)
library(ggplot2)
library(cowplot)



# in/out
io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data"


io$epi_ect <- "/bi/home/clarks/gastrulation/plots/motifs/diffexp_diffacc/background_all_k27ac/260419/down/E4.5EPI_vs_E7.5Ect_H3K27ac_distal_E7.5_Ect_intersect12.tsv"
io$ect_mesend <- "/bi/home/clarks/gastrulation/plots/motifs/diffexp_diffacc/background_all_k27ac/260419/up/E7.5Ect_vs_E7.5EndMes_H3K27ac_distal_E7.5_Ect_intersect12.tsv"



toplot <- map(list(io$epi_ect, io$ect_mesend), fread) %>%
  map(dcast, tf + motif + lab ~ comparison, value.var = c("log_fdr", "logFC")) %>%
  map(setkey, tf, motif, lab) %>%
  purrr::reduce(merge)


ggplot(toplot, aes(log_fdr_E4.5EPI_vs_E7.5Ect, log_fdr_E7.5Ect_vs_E7.5EndMes, colour = logFC_E4.5EPI_vs_E7.5Ect)) +
  geom_point()
