library(data.table)
library(purrr)
library(ggplot2)
library(cowplot)
library(ggpubr)


io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data/"
io$motif_acc <- "/bi/scratch/Stephen_Clark/gastrulation_data/acc/parsed/fimo_motifs/020519/200bp/"
io$meta <- "/bi/scratch/Stephen_Clark/gastrulation_data/sample_metadata.txt"


opts <- list()
opts$tfs <- c("POU5F1", "KLF4", "SP8", "Sox2")
opts$lineages <- c("Epiblast", "Primitive_endoderm", "Ectoderm", "Mesoderm")
# opts$lineages <- c("Epiblast", "Ectoderm")
opts$anno <- "E7.5_Ect_intersect12"

# fun
fread_gz <- function(filename, ...){
  f <- file(filename)
  type <- summary(f)$class
  close.connection(f)
  if (type == "gzfile") {
    filename <- paste("zcat", filename)
    return(fread(cmd = filename, ...))
  }
  fread(filename, ...)
}

fwrite_tsv <- partial(fwrite, sep = "\t", na = "NA")


################################################################################


meta <- fread(io$meta) %>%
  .[pass_accQC == TRUE] %>%
  .[, lineage := lineage10x_2] %>%
  .[lineage %in% opts$lineages]

acc <- paste0(io$motif_acc, "/", opts$tfs, ".tsv.gz") %>%
  .[file.exists(.)] %>%
  map(fread_gz) %>%
  rbindlist() %>%
  .[id_motif %like% opts$anno] %>%
  merge(meta[, .(sample, stage_lineage = paste0(stage, "_", lineage))], by = "sample")

toplot_cells <- acc[, .(rate = mean(rate)), .(tf, sample, stage_lineage)]
toplot_loci  <- acc[, .(rate = mean(rate)), .(tf, id_motif, stage_lineage)]

ggplot(toplot_cells, aes(stage_lineage, rate, fill = stage_lineage)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~tf) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position = "none")



